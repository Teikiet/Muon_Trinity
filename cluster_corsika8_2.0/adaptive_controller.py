#!/usr/bin/env python3
"""
Adaptive controller for Trinity muon simulations.

State machine: shell progression -> zenith branches -> height cells -> N-doubling.

Each invocation does exactly one cycle:
  1. Load/init state (flock-protected)
  2. Ingest per-zenith CSVs from previous cycle's pending jobs
  3. Advance state machine (normalize branches, maybe advance shell)
  4. Plan next work items (cells needing more seeds)
  5. Respect --max-queued
  6. sbatch sim jobs
  7. sbatch one save job (dependency on all sims)
  8. sbatch next controller (dependency on save job)
  9. Record pending, write state, exit
"""

import argparse
import csv
import fcntl
import glob
import json
import os
import shlex
import subprocess
import sys
import tempfile
import time
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
SAVE_DIR = SCRIPT_DIR / "save_data2csv"
if str(SAVE_DIR) not in sys.path:
    sys.path.insert(0, str(SAVE_DIR))

from save_CARE2csv_chunk_tree import (  # noqa: E402
    _fmt,
    _fmt_angle,
    _fmt_offset,
    make_correction_paths,
    make_path,
    make_run_dir,
)


# ── Constants ──────────────────────────────────────────────────────────────
BASE_PATH = Path("/scratch/general/vast/u1520754/muon_sim_chain_tree_adaptive")
CONDA_SH = "/uufs/chpc.utah.edu/common/home/u1520754/miniconda3/etc/profile.d/conda.sh"
WORKER_SCRIPT = SAVE_DIR / "save_CARE2csv_chunk_tree.py"
APPEND_SCRIPT = SAVE_DIR / "append_adaptive_rows.py"
SLURM_SCRIPT = SCRIPT_DIR / "base_muon_event_generation" / "run_corsika8_trinity_chain_tree.slurm"
DEFAULT_HADRON_MODEL = "SIBYLL-2.3d"
DEFAULT_OBS_LEVEL = 2944
MAX_SUBMIT_RETRIES = 20
RETRY_SLEEP_SECONDS = 30


# ── Utility helpers ────────────────────────────────────────────────────────
def compute_cherenkov_radius(_energy_string):
    """Hard-coded Cherenkov radius; mirrors submit_multiple_continuous_tree.sh."""
    return 15.0


def normalize_seed(value):
    return int(float(value))


def format_zen(value):
    return f"{float(value):.1f}"


def format_height(value):
    return str(int(round(float(value))))


def _offset_key(value):
    """Match _fmt_offset from save_CARE2csv_chunk_tree.py."""
    return _fmt_offset(value)


def energy_output_dir(pid, energy_string, radius):
    return BASE_PATH / f"Muon_pid{pid}_E{energy_string}_R{radius}"


def csv_dir_for(pid, energy_string, radius):
    return energy_output_dir(pid, energy_string, radius) / "csv_output"


def state_path_for(pid, energy_string, radius):
    return energy_output_dir(pid, energy_string, radius) / "adaptive_state.json"


# ── Shell / branch helpers ──────────────────────────────────────────────────
def shell_branch_keys(shell_index, config):
    """Return the zenith key(s) for a given shell, respecting zen_max."""
    zen_start = float(config["zen_start"])
    delta_zen = float(config["delta_zen"])
    zen_max = float(config.get("zen_max", 92.0))
    if shell_index == 0:
        return [format_zen(zen_start)]
    low = zen_start - shell_index * delta_zen
    high = zen_start + shell_index * delta_zen
    keys = []
    if low <= zen_max:
        keys.append(format_zen(low))
    if high <= zen_max:
        keys.append(format_zen(high))
    return keys


def current_shell_zeniths(state):
    return shell_branch_keys(int(state["current_shell"]), state["config"])


def current_shell_branches(state):
    """Yield (branch_key, branch) for the active shell."""
    for k in current_shell_zeniths(state):
        branch = state["branches"].get(k)
        if branch is not None:
            yield k, branch


def ensure_shell_structure(state, shell_index):
    """Create branches for this shell if missing."""
    keys = shell_branch_keys(shell_index, state["config"])
    for k in keys:
        if k not in state["branches"]:
            state["branches"][k] = {
                "status": "active",
                "next_h": int(state["config"]["h_start"]),
                "cells": {},
                "trigger_count": 0,
                "shell": shell_index,
            }


def current_group_hs(branch, config):
    """Return the next batch of h values for this branch."""
    h_values = config.get("h_values") or []
    if h_values:
        # Custom height list mode: advance by index into the list.
        next_idx = int(branch.get("next_h_idx", 0))
        batch = int(config["h_batch"])
        hs = []
        for i in range(batch):
            idx = next_idx + i
            if idx >= len(h_values):
                break
            hs.append(h_values[idx])
        return hs
    # Arithmetic mode (backward compatible)
    next_h = int(branch.get("next_h", config["h_start"]))
    batch = int(config["h_batch"])
    delta = int(config["delta_h"])
    h_max = int(config["h_max"])
    hs = []
    for i in range(batch):
        h = next_h + i * delta
        if h > h_max:
            break
        hs.append(h)
    return hs


def parse_h_values(value):
    """Parse a comma-separated string of heights into a sorted, deduplicated list of ints."""
    if not value:
        return []
    out = []
    seen = set()
    for part in str(value).split(","):
        part = part.strip()
        if not part:
            continue
        try:
            h = int(round(float(part)))
        except (ValueError, TypeError):
            continue
        if h not in seen:
            seen.add(h)
            out.append(h)
    return sorted(out)


def default_config(args):
    return {
        "pid": args.pid,
        "radius": args.radius,
        "tel_y": args.tel_y,
        "n_init": args.n_init,
        "n_max": args.n_max,
        "min_triggers": getattr(args, "min_triggers", 1),
        "h_start": args.h_start,
        "h_max": args.h_max,
        "delta_h": args.delta_h,
        "h_terminate_above": args.h_terminate_above,
        "h_batch": args.h_batch,
        "zen_start": args.zen_start,
        "zen_max": getattr(args, "zen_max", 92.0),
        "delta_zen": args.delta_zen,
        "az": args.az,
        "tel_x": args.tel_x,
        "tel_z": args.tel_z,
        "trigger_pe": args.trigger_pe,
        "max_queued": args.max_queued,
        "save_batches": getattr(args, "save_batches", 1),
        "obs_level": getattr(args, "obs_level", DEFAULT_OBS_LEVEL),
        "hadron_model": getattr(args, "hadron_model", DEFAULT_HADRON_MODEL),
        "update_time_stamp": getattr(args, "update_time_stamp", False),
        "delete_log_on_success": getattr(args, "delete_log_on_success", False),
        "energy_string": args.energy,
        "h_values": parse_h_values(getattr(args, "h_values", None)),
    }


# ── State file I/O ────────────────────────────────────────────────────────
def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def write_json_atomic(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", delete=False, dir=str(path.parent), encoding="utf-8") as f:
        json.dump(data, f, indent=2)
        f.write("\n")
        tmp = Path(f.name)
    os.replace(tmp, path)


def acquire_state_lock(state_path):
    lock_path = Path(f"{state_path}.lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    h = open(lock_path, "a+")
    fcntl.flock(h, fcntl.LOCK_EX)
    return h


def release_state_lock(h):
    try:
        fcntl.flock(h, fcntl.LOCK_UN)
    finally:
        h.close()


# ── CSV / cell helpers ─────────────────────────────────────────────────────
def load_csv_rows(csv_path):
    """Load rows from a CSV file; return empty list if missing."""
    path = Path(csv_path)
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return list(reader)


def cell_from_rows(rows, trigger_pe):
    """
    From a list of CSV rows for one cell, determine:
      (n_success, n_failed, trigger_count)
    """
    n_success = 0
    n_failed = 0
    trigger_count = 0
    for row in rows:
        file_found = str(row.get("file_found", "0")).strip()
        if file_found == "1":
            n_success += 1
            try:
                max_pe = float(row.get("max_pe", 0))
                if max_pe >= float(trigger_pe):
                    trigger_count += 1
            except (ValueError, TypeError):
                pass
        else:
            n_failed += 1
    return n_success, n_failed, trigger_count


def choose_target(n_success, n_init, n_max):
    """
    Return the next N-doubling target, or None if already at n_max.
    Cumulative doubling: n_init, 2*n_init, 4*n_init, ...
    """
    target = n_init
    while target <= n_max:
        if n_success < target:
            return target
        target *= 2
    return None


# ── CSV update ─────────────────────────────────────────────────────────────
def update_state_from_csvs(state):
    """
    Scan per-zenith CSVs for the current energy, update cell counters.
    Returns number of rows seen.
    """
    config = state["config"]
    pid = state["pid"]
    energy = state.get("energy_string", config.get("energy_string", ""))
    radius = state["radius"]
    tel_y = state["tel_y"]
    csv_root = csv_dir_for(pid, energy, radius)
    pattern = str(csv_root / f"adaptive_pid{pid}_E{energy}_R{radius}_y{tel_y}_zen*.csv")

    rows_seen = 0
    for csv_path in sorted(Path(p) for p in glob.glob(pattern)):
        rows = load_csv_rows(csv_path)
        if not rows:
            continue
        rows_seen += len(rows)

        # Extract zenith from filename: adaptive_pid13_E1.12202e1_R5_y0_zen90.0.csv
        fname = csv_path.name
        parts = fname.rsplit("_zen", 1)
        if len(parts) != 2:
            continue
        zen_key = parts[1].replace(".csv", "")
        branch = state["branches"].setdefault(
            zen_key,
            {
                "status": "active",
                "next_h": int(config["h_start"]),
                "cells": {},
                "trigger_count": 0,
                "shell": int(state.get("current_shell", 0)),
            },
        )
        cells = branch.setdefault("cells", {})

        # Group rows by height
        cell_groups = {}
        for row in rows:
            ck = format_height(row.get("height", 0))
            cell_groups.setdefault(ck, []).append(row)

        for cell_key, cell_rows in cell_groups.items():
            cell = cells.setdefault(
                cell_key,
                {
                    "n_success": 0,
                    "n_failed": 0,
                    "triggered": False,
                    "status": "active",
                    "attempts": {},
                },
            )
            # Per-cell fields (Change 3) — optional on read so an existing
            # state file without them still loads (defaults below).
            cell.setdefault("target", int(config["n_init"]))
            cell.setdefault("extra_indices", [])
            cell.setdefault("seeds", [])

            # Collect every attempted seed index for this cell (all rows,
            # regardless of file_found). This is the `already` set used to
            # derive per-cell cumulative indices (Change 2).
            attempted = set()
            for row in cell_rows:
                try:
                    attempted.add(int(float(row.get("seed", 0))))
                except (ValueError, TypeError):
                    pass
            cell["seeds"] = sorted(attempted)

            n_success, n_failed, trigger_count = cell_from_rows(cell_rows, config["trigger_pe"])
            cell["n_success"] = n_success
            cell["n_failed"] = n_failed
            cell["trigger_count"] = trigger_count
            min_triggers = int(config.get("min_triggers", 1))
            if trigger_count >= min_triggers:
                cell["triggered"] = True
                cell["status"] = "triggered"
            elif cell.get("status") not in {"triggered", "failed"}:
                if n_success >= int(config["n_max"]):
                    cell["status"] = "exhausted"
                elif cell.get("status") not in {"triggered", "exhausted"}:
                    cell["status"] = "active"

        branch["trigger_count"] = sum(
            1 for c in cells.values() if c.get("triggered")
        )

    return rows_seen


# ── State machine ──────────────────────────────────────────────────────────
def branch_has_trigger(branch):
    """Shell termination check: any cell with >= 1 trigger (regardless of min_triggers)."""
    return any(c.get("trigger_count", 0) >= 1 for c in branch.get("cells", {}).values())


def branch_group_terminal(branch, hs, h_terminate_above):
    """Check if all h-values in the current group are terminal."""
    cells = branch.get("cells", {})
    if not hs:
        return True, False
    terminal = True
    exhausted_high = False
    for h_value in hs:
        ck = format_height(h_value)
        cell = cells.get(ck)
        if cell is None:
            return False, False
        status = cell.get("status", "active")
        if status == "active":
            terminal = False
        if (
            status == "exhausted"
            and not cell.get("triggered")
            and float(h_value) >= float(h_terminate_above)
        ):
            exhausted_high = True
    return terminal, exhausted_high


def normalize_active_branches(state):
    """
    Advance branches past cells that are already complete/terminal.
    """
    config = state["config"]
    h_terminate_above = float(config.get("h_terminate_above", 20000))
    h_values = config.get("h_values") or []
    for branch in state["branches"].values():
        if branch.get("status") != "active":
            continue
        while branch.get("status") == "active":
            hs = current_group_hs(branch, config)
            if not hs:
                branch["status"] = "complete"
                break
            terminal, ex_high = branch_group_terminal(branch, hs, h_terminate_above)
            if not terminal:
                break
            if ex_high:
                branch["status"] = "terminated_early"
                break
            if h_values:
                next_idx = int(branch.get("next_h_idx", 0)) + int(config["h_batch"])
                if next_idx >= len(h_values):
                    branch["status"] = "complete"
                    break
                branch["next_h_idx"] = next_idx
            else:
                next_h = int(branch["next_h"]) + int(config["h_batch"]) * int(config["delta_h"])
                if next_h > int(config["h_max"]):
                    branch["status"] = "complete"
                    break
                branch["next_h"] = next_h


def current_shell_finished(state):
    keys = current_shell_zeniths(state)
    branches = [state["branches"].get(k) for k in keys]
    branches = [b for b in branches if b is not None]
    if not branches:
        return False
    return all(b.get("status") != "active" for b in branches)


def current_shell_any_trigger(state):
    for _, branch in current_shell_branches(state):
        if branch_has_trigger(branch):
            return True
    return False


def maybe_advance_shell(state):
    """If current shell finished, check triggers and advance or done."""
    if not current_shell_finished(state):
        return False
    shell_index = int(state["current_shell"])
    any_trigger = current_shell_any_trigger(state)
    state.setdefault("shell_history", []).append({
        "shell": shell_index,
        "zeniths": current_shell_zeniths(state),
        "any_trigger": bool(any_trigger),
    })
    if not any_trigger:
        state["status"] = "done"
        return True
    next_shell = shell_index + 1
    next_keys = shell_branch_keys(next_shell, state["config"])
    if not next_keys:
        # All branches of the next shell would exceed zen_max — stop expanding.
        state["status"] = "done"
        return True
    state["current_shell"] = next_shell
    ensure_shell_structure(state, next_shell)
    return True


def apply_attempt_rounds(state, submitted_jobs):
    """Track one attempt round per (cell, target); cap at 5 rounds per target."""
    seen = set()
    for job in submitted_jobs:
        branch_key = job.get("branch_key", format_zen(job["zen"]))
        height_key = format_height(job["h"])
        target = int(job.get("target", 0))
        key = (branch_key, height_key, target)
        if key in seen:
            continue
        seen.add(key)
        branch = state["branches"].get(branch_key)
        if branch is None:
            continue
        cell = branch.get("cells", {}).get(height_key)
        if cell is None:
            continue
        attempts = cell.setdefault("attempts", {})
        attempts[str(target)] = int(attempts.get(str(target), 0)) + 1


# ── Job planning ───────────────────────────────────────────────────────────
def plan_jobs_for_state(state):
    """
    Determine which cells need more seeds.

    Seed indices are per-cell cumulative: a cell always draws from the
    contiguous range 1..target, where target is its current N-doubling
    target. `already` = every seed index already attempted for the cell
    (from the per-zenith CSV, all rows regardless of file_found), so a
    failed index is never re-derived within the range. If the fixed 1..target
    set is exhausted but the cell is still short of `target` successful runs
    (permanently-broken index), the range is extended upward via
    `extra_indices` (Change 2), keeping the N-successful-runs contract.

    Returns (jobs, branch_plan, cell_plan) where jobs is a list of dicts,
    branch_plan records (branch_key, hs), and cell_plan records per-cell
    (zen, h, target, already, needed) for dry-run reporting.
    """
    config = state["config"]
    n_init = int(config["n_init"])
    n_max = int(config["n_max"])
    jobs = []
    branch_plan = []
    cell_plan = []

    for branch_key, branch in current_shell_branches(state):
        if branch.get("status") != "active":
            continue
        hs = current_group_hs(branch, config)
        if not hs:
            continue
        branch_plan.append((branch_key, list(hs)))

        for h_value in hs:
            ck = format_height(h_value)
            cell = branch["cells"].setdefault(
                ck,
                {
                    "n_success": 0,
                    "n_failed": 0,
                    "triggered": False,
                    "status": "active",
                    "attempts": {},
                    "target": n_init,
                    "extra_indices": [],
                    "seeds": [],
                },
            )
            if cell.get("triggered") or cell.get("status") in ("failed", "exhausted"):
                continue
            target = choose_target(cell.get("n_success", 0), n_init, n_max)
            if target is None:
                cell["status"] = "exhausted"
                continue
            cell["target"] = target

            attempts = int(cell.get("attempts", {}).get(str(target), 0))
            already = set(int(s) for s in (cell.get("seeds", []) or []))
            extra_indices = [int(s) for s in (cell.get("extra_indices", []) or [])]
            extra_set = set(extra_indices)

            # Untried indices within the fixed 1..target range.
            needed = [k for k in range(1, target + 1)
                      if k not in already and k not in extra_set]
            n_success = int(cell.get("n_success", 0))
            shortfall = max(0, target - n_success)

            if shortfall <= 0:
                continue

            # Fixed 1..target set fully consumed but still short of target
            # successes (some indices permanently failed). Extend upward with
            # extra indices until len(needed) + n_success == target. This also
            # fires once the 5-round replacement cap is reached while short.
            fixed_exhausted = (not needed) and all(
                k in already or k in extra_set for k in range(1, target + 1)
            )
            if (attempts >= 5 or fixed_exhausted) and len(needed) + n_success < target:
                n_extra = target - n_success - len(needed)
                nxt = (max(extra_indices) if extra_indices else target) + 1
                added = []
                while n_extra > 0:
                    if nxt not in already and nxt not in extra_set:
                        added.append(nxt)
                        extra_set.add(nxt)
                        n_extra -= 1
                    nxt += 1
                if added:
                    cell["extra_indices"] = sorted(extra_indices + added)
                    needed.extend(added)

            if attempts >= 5:
                # Retain the 5-round cap and the status="failed" outcome.
                cell["status"] = "failed"

            cell_plan.append({
                "zen": float(branch_key),
                "h": int(h_value),
                "target": target,
                "already": sorted(already),
                "needed": list(needed),
            })

            if not needed:
                continue
            for k in needed:
                jobs.append({
                    "branch_key": branch_key,
                    "zen": float(branch_key),
                    "az": float(config["az"]),
                    "h": int(h_value),
                    "tel_x": float(config["tel_x"]),
                    "tel_z": float(config["tel_z"]),
                    "seed": k,
                    "target": target,
                    "height_key": ck,
                })
    return jobs, branch_plan, cell_plan


# ── Build command helpers ──────────────────────────────────────────────────
def build_sim_command(
    output_base_dir, pid, energy_string, radius, tel_y,
    zenith, azimuth, height, tel_x, tel_z, seed, config,
    rerun_mode=False,
):
    """Build the sbatch command for one simulation job.

    Mirrors the convention in submit_multiple_continuous_tree.sh exactly.
    """
    log_dir = Path(output_base_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"s{seed}_zen{format_zen(zenith)}_az{format_zen(azimuth)}_h{format_height(height)}_x{_fmt_offset(tel_x)}_z{_fmt_offset(tel_z)}.log"
    cherenkov_radius = compute_cherenkov_radius(energy_string)
    rerun_event = "detector_sim" if rerun_mode else "false"

    command = [
        "sbatch",
        "--parsable",
        f"--job-name=corsika8_trinity_pid{pid}_E{energy_string}_s{seed}",
        str(SLURM_SCRIPT),
        str(output_base_dir),
        str(log_file),
        str(pid),
        str(energy_string),
        format_zen(zenith),
        format_zen(azimuth),
        format_height(height),
        str(config["obs_level"]),
        str(tel_x),
        str(tel_y),
        str(tel_z),
        f"{cherenkov_radius:.1f}",
        str(radius),
        str(seed),
        str(config["hadron_model"]),
        rerun_event,
        f"--update_time_stamp={str(config['update_time_stamp']).lower()}",
    ]
    if config.get("delete_log_on_success"):
        command.append("--delete_log_on_success")
    return command


def build_save_command(chunk_file, result_csv, config, mode="append"):
    """Build the sbatch command for the save job.

    Runs save_CARE2csv_chunk_tree.py then append_adaptive_rows.py.
    """
    save_log_dir = Path.home() / "csv_logs" / f"adaptive_E{config['energy_string']}"
    save_log_dir.mkdir(parents=True, exist_ok=True)

    wrap = (
        f"source {shlex.quote(CONDA_SH)} && "
        f"conda activate jupyter_env && "
        f"python3 {shlex.quote(str(WORKER_SCRIPT))} "
        f"--chunk-file {shlex.quote(str(chunk_file))} "
        f"--output {shlex.quote(str(result_csv))} "
        f"--pid {config['pid']} --energy-str {shlex.quote(str(config['energy_string']))} "
        f"--radius {config['radius']} --tel-y {config['tel_y']} "
        f"--base-path {shlex.quote(str(BASE_PATH))} "
        f"--correction-report-name metadata.yaml "
        f"--Max_PE_cut {config['trigger_pe']} && "
        f"python3 {shlex.quote(str(APPEND_SCRIPT))} "
        f"--result-csv {shlex.quote(str(result_csv))} "
        f"--pid {config['pid']} --energy-str {shlex.quote(str(config['energy_string']))} "
        f"--radius {config['radius']} --tel-y {config['tel_y']} "
        f"--base-path {shlex.quote(str(BASE_PATH))} --mode {mode}"
    )

    return [
        "sbatch",
        "--parsable",
        f"--job-name=adaptive_save_E{config['energy_string']}",
        f"--output={save_log_dir / 'save_%j.out'}",
        f"--error={save_log_dir / 'save_%j.err'}",
        "--account=owner-guest",
        "--partition=kingspeak-guest",
        "--time=1:30:00",
        "--mem=4G",
        "--wrap",
        wrap,
    ]


def build_save_worker_command(chunk_file, result_csv, config):
    """Build the sbatch command for ONE batched save worker job.

    Processes a chunk of combos and writes a temp result CSV (no append).
    """
    save_log_dir = Path.home() / "csv_logs" / f"adaptive_E{config['energy_string']}"
    save_log_dir.mkdir(parents=True, exist_ok=True)

    wrap = (
        f"source {shlex.quote(CONDA_SH)} && "
        f"conda activate jupyter_env && "
        f"python3 {shlex.quote(str(WORKER_SCRIPT))} "
        f"--chunk-file {shlex.quote(str(chunk_file))} "
        f"--output {shlex.quote(str(result_csv))} "
        f"--pid {config['pid']} --energy-str {shlex.quote(str(config['energy_string']))} "
        f"--radius {config['radius']} --tel-y {config['tel_y']} "
        f"--base-path {shlex.quote(str(BASE_PATH))} "
        f"--correction-report-name metadata.yaml "
        f"--Max_PE_cut {config['trigger_pe']}"
    )

    return [
        "sbatch",
        "--parsable",
        f"--job-name=adaptive_save_E{config['energy_string']}",
        f"--output={save_log_dir / 'save_%j.out'}",
        f"--error={save_log_dir / 'save_%j.err'}",
        "--account=owner-guest",
        "--partition=kingspeak-guest",
        "--time=1:30:00",
        "--mem=4G",
        "--wrap",
        wrap,
    ]


def build_merge_command(batch_result_csvs, config, mode="append"):
    """Build the sbatch command that merges several batch result CSVs into the
    per-zenith destination CSVs in a single flocked append/replace."""
    save_log_dir = Path.home() / "csv_logs" / f"adaptive_E{config['energy_string']}"
    save_log_dir.mkdir(parents=True, exist_ok=True)

    rc_args = " ".join(
        f"--result-csv {shlex.quote(str(p))}" for p in batch_result_csvs
    )
    wrap = (
        f"source {shlex.quote(CONDA_SH)} && "
        f"conda activate jupyter_env && "
        f"python3 {shlex.quote(str(APPEND_SCRIPT))} {rc_args} "
        f"--pid {config['pid']} --energy-str {shlex.quote(str(config['energy_string']))} "
        f"--radius {config['radius']} --tel-y {config['tel_y']} "
        f"--base-path {shlex.quote(str(BASE_PATH))} --mode {mode}"
    )

    return [
        "sbatch",
        "--parsable",
        f"--job-name=adaptive_merge_E{config['energy_string']}",
        f"--output={save_log_dir / 'merge_%j.out'}",
        f"--error={save_log_dir / 'merge_%j.err'}",
        "--account=owner-guest",
        "--partition=kingspeak-guest",
        "--time=0:30:00",
        "--mem=2G",
        "--wrap",
        wrap,
    ]


def build_next_controller_command(args, state_file, save_job_id=None, begin_delay=False):
    """Build the sbatch command for the next controller cycle."""
    log_dir = Path.home() / "csv_logs" / f"adaptive_E{args.energy}"
    log_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "sbatch",
        "--parsable",
        f"--job-name=adaptive_controller_E{args.energy}",
        f"--output={log_dir / 'cycle_%j.out'}",
        f"--error={log_dir / 'cycle_%j.err'}",
        "--account=owner-guest",
        "--partition=kingspeak-guest",
        "--time=0:15:00",
        "--mem=2G",
    ]
    if begin_delay:
        cmd.append("--begin=now+2minutes")
    if save_job_id:
        cmd.append(f"--dependency=afterany:{save_job_id}")

    # Build the python command for --wrap
    py_args = [str(SCRIPT_DIR / "adaptive_controller.py")]
    for attr in [
        "energy", "n_init", "n_max", "min_triggers", "h_start", "h_max", "delta_h",
        "h_terminate_above", "h_batch", "zen_start", "zen_max", "delta_zen", "az",
        "trigger_pe", "max_queued", "save_batches", "pid", "radius", "tel_y",
        "obs_level", "hadron_model", "h_values",
    ]:
        val = getattr(args, attr, None)
        if val is not None:
            py_args.append(f"--{attr.replace('_', '-')}={val}")
    if args.update_time_stamp:
        py_args.append("--update-time-stamp")
    if args.delete_log_on_success:
        py_args.append("--delete-log-on-success")
    if args.resume:
        py_args.append("--resume")
    py_args.append(f"--state-file={state_file}")

    wrap = (
        f"source {shlex.quote(CONDA_SH)} && "
        f"conda activate jupyter_env && "
        f"python3 {' '.join(shlex.quote(a) for a in py_args)}"
    )
    cmd.extend(["--wrap", wrap])
    return cmd


def build_chunk_file(work_dir, jobs, name="cycle_chunk.json"):
    """Write the chunk JSON file for this cycle's jobs (or one batch of them)."""
    chunk_file = Path(work_dir) / name
    chunk_file.parent.mkdir(parents=True, exist_ok=True)
    payload = [
        [job["zen"], job["az"], job["h"], job["tel_x"], job["tel_z"], job["seed"]]
        for job in jobs
    ]
    with chunk_file.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
        f.write("\n")
    return chunk_file


def split_list(items, n):
    """Split items into n contiguous (non-empty where possible) slices."""
    if n <= 0 or not items:
        return [list(items)] if n > 0 else []
    base = len(items) // n
    rem = len(items) % n
    out = []
    i = 0
    for j in range(n):
        size = base + (1 if j < rem else 0)
        out.append(items[i:i + size])
        i += size
    return out


# ── Submit with retry ──────────────────────────────────────────────────────
def submit_with_retry(command, dry_run=False):
    """Submit a job via sbatch, retrying on QOSMaxSubmitJobPerUserLimit."""
    if dry_run:
        print("DRY-RUN:", shlex.join(command))
        return None

    last_error = ""
    for attempt in range(1, MAX_SUBMIT_RETRIES + 1):
        proc = subprocess.run(command, capture_output=True, text=True)
        output = (proc.stdout or "").strip()
        error = (proc.stderr or "").strip()
        last_error = "\n".join(p for p in [output, error] if p)
        if proc.returncode == 0 and output:
            return output.split()[-1]
        if (
            "QOSMaxSubmitJobPerUserLimit" in last_error
            or "Resource temporarily unavailable" in last_error
            or "Socket timed out" in last_error
        ):
            print(f"  [retry {attempt}/{MAX_SUBMIT_RETRIES}] sbatch limit hit, sleeping {RETRY_SLEEP_SECONDS}s")
            time.sleep(RETRY_SLEEP_SECONDS)
            continue
        raise RuntimeError(f"sbatch failed: {last_error}")
    raise RuntimeError(f"sbatch failed after {MAX_SUBMIT_RETRIES} retries: {last_error}")


# ── Rerun failed mode ──────────────────────────────────────────────────────
def simulate_rerun_failed(args):
    """Re-submit failed (file_found=0) rows with the same seeds."""
    config = default_config(args)
    config["energy_string"] = args.energy
    csv_root = csv_dir_for(args.pid, args.energy, args.radius)
    pattern = str(csv_root / f"adaptive_pid{args.pid}_E{args.energy}_R{args.radius}_y{args.tel_y}_zen*.csv")
    failed = []
    seen = set()
    for csv_path in sorted(Path(p) for p in glob.glob(pattern)):
        for row in load_csv_rows(csv_path):
            if str(row.get("file_found", "0")).strip() != "0":
                continue
            key = (
                format_zen(row.get("zen", 0.0)),
                format_zen(row.get("az", 0.0)),
                format_height(row.get("height", 0)),
                _fmt_offset(row.get("tel_x", 0.0)),
                _fmt_offset(row.get("tel_z", 0.0)),
                normalize_seed(row.get("seed", 0)),
            )
            if key in seen:
                continue
            seen.add(key)
            failed.append({
                "branch_key": key[0],
                "zen": float(key[0]),
                "az": float(key[1]),
                "h": int(key[2]),
                "tel_x": float(key[3]),
                "tel_z": float(key[4]),
                "seed": int(key[5]),
                "target": 0,
            })
    if not failed:
        print("No failed rows found to rerun.")
        return 0

    work_dir = energy_output_dir(args.pid, args.energy, args.radius) / "adaptive_work" / "rerun_failed"
    work_dir.mkdir(parents=True, exist_ok=True)
    chunk_file = build_chunk_file(work_dir, failed)
    result_csv = work_dir / "rerun_failed_result.csv"

    sim_jids = []
    for job in failed:
        cmd = build_sim_command(
            energy_output_dir(args.pid, args.energy, args.radius),
            args.pid, args.energy, args.radius, args.tel_y,
            job["zen"], job["az"], job["h"], job["tel_x"], job["tel_z"], job["seed"],
            config, rerun_mode=True,
        )
        jid = submit_with_retry(cmd, dry_run=args.dry_run)
        if jid is not None:
            sim_jids.append(jid)

    if args.dry_run:
        save_cmd = build_save_command(chunk_file, result_csv, config, mode="replace")
        print("DRY-RUN:", shlex.join(save_cmd))
        print(f"Would retry {len(failed)} failed rows.")
        return 0

    if sim_jids:
        save_cmd = build_save_command(chunk_file, result_csv, config, mode="replace")
        wrap_idx = save_cmd.index("--wrap")
        save_cmd.insert(wrap_idx, f"--dependency=afterany:{':'.join(sim_jids)}")
        save_jid = submit_with_retry(save_cmd)
        print(f"Retried {len(failed)} failed rows; save job {save_jid}.")
    else:
        print(f"Retried {len(failed)} failed rows; no sims submitted.")
    return 0


# ── Build initial state ───────────────────────────────────────────────────
def build_state_template(args):
    config = default_config(args)
    state = {
        "version": 2,
        "energy_string": args.energy,
        "pid": args.pid,
        "radius": args.radius,
        "tel_y": args.tel_y,
        "config": config,
        "status": "active",
        "cycle": 0,
        "current_shell": 0,
        "branches": {},
        "shell_history": [],
        "pending": [],
        "last_save_job": None,
    }
    ensure_shell_structure(state, 0)
    return state


def guard_old_schema(state):
    """Reject a state file that still uses the removed global-seed schema."""
    if "seed_counter" in state:
        print(
            "ERROR: state file uses the old global-seed schema (seed_counter present).\n"
            "Per-cell seed indices are incompatible with it. Delete the output tree and\n"
            "state file, then relaunch. Aborting.",
            file=sys.stderr,
        )
        raise SystemExit(1)


# ── Main ───────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Adaptive controller for Trinity muon simulations.")
    parser.add_argument("--energy", required=True)
    parser.add_argument("--state-file")
    parser.add_argument("--pid", type=int, default=13)
    parser.add_argument("--radius", type=int, default=5)
    parser.add_argument("--tel-y", type=float, default=0)
    parser.add_argument("--obs-level", type=int, default=DEFAULT_OBS_LEVEL)
    parser.add_argument("--hadron-model", default=DEFAULT_HADRON_MODEL)
    parser.add_argument("--n-init", type=int, default=1)
    parser.add_argument("--n-max", type=int, default=64)
    parser.add_argument("--min-triggers", type=int, default=1,
                        help="Minimum number of triggered events (max_pe >= trigger_pe) "
                             "required before a cell stops N-doubling (default 1). "
                             "Shell termination still uses >= 1 trigger.")
    parser.add_argument("--h-start", type=int, default=3000)
    parser.add_argument("--h-max", type=int, default=50000)
    parser.add_argument("--delta-h", type=int, default=200)
    parser.add_argument("--h-terminate-above", type=int, default=20000)
    parser.add_argument("--h-batch", type=int, default=1)
    parser.add_argument("--zen-start", type=float, default=90.0)
    parser.add_argument("--zen-max", type=float, default=92.0,
                        help="Maximum zenith (deg). The scan stops expanding "
                             "when a shell's branches would exceed this value.")
    parser.add_argument("--h-values", default=None,
                        help="Comma-separated list of injection heights (m). "
                             "When set, replaces h_start/delta_h/h_max scanning "
                             "with this explicit list.")
    parser.add_argument("--delta-zen", type=float, default=0.3)
    parser.add_argument("--az", type=float, default=270.0)
    parser.add_argument("--tel-x", type=float, default=0.0)
    parser.add_argument("--tel-z", type=float, default=0.0)
    parser.add_argument("--trigger-pe", type=float, default=20.0)
    parser.add_argument("--max-queued", type=int, default=900)
    parser.add_argument("--save-batches", type=int, default=1,
                        help="Split each cycle's save into N parallel worker jobs "
                             "followed by one merge job (default 1 = single save).")
    parser.add_argument("--update-time-stamp", action="store_true")
    parser.add_argument("--delete-log-on-success", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--rerun-failed", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    if args.rerun_failed:
        return simulate_rerun_failed(args)

    state_file = Path(args.state_file) if args.state_file else state_path_for(
        args.pid, args.energy, args.radius,
    )

    # ── Dry-run: show planned work without submitting ──
    if args.dry_run:
        if state_file.exists():
            state = load_json(state_file)
        else:
            state = build_state_template(args)
        guard_old_schema(state)
        state["config"] = default_config(args)
        update_state_from_csvs(state)
        normalize_active_branches(state)
        maybe_advance_shell(state)
        if state.get("status") == "done":
            print("DRY-RUN: scan already complete.")
            return 0
        jobs, bp, cell_plan = plan_jobs_for_state(state)
        print(f"DRY-RUN: current shell={state['current_shell']}, zeniths={current_shell_zeniths(state)}")
        print(f"DRY-RUN: would submit {len(jobs)} sim jobs")
        # Per-cell seed index plan (Change 7)
        for cp in cell_plan:
            already_repr = "{" + ",".join(str(s) for s in cp["already"]) + "}"
            print(f"CELL zen={cp['zen']:.1f} h={cp['h']}  target={cp['target']}  "
                  f"already={already_repr} needed={cp['needed']}")
        if jobs:
            print(f"DRY-RUN: first job = zen={jobs[0]['zen']:.1f} h={jobs[0]['h']} seed={jobs[0]['seed']}")
            # Print the literal sim sbatch command for the first job
            output_base = energy_output_dir(args.pid, args.energy, args.radius)
            sim_cmd = build_sim_command(
                output_base, args.pid, args.energy, args.radius, args.tel_y,
                jobs[0]["zen"], jobs[0]["az"], jobs[0]["h"],
                jobs[0]["tel_x"], jobs[0]["tel_z"], jobs[0]["seed"],
                state["config"],
            )
            print(f"DRY-RUN: literal sim sbatch command:")
            print(f"DRY-RUN: {shlex.join(sim_cmd)}")
            # List positional arguments in order with cycle-0 values
            j0 = jobs[0]
            log_file_str = str(Path(output_base) / "logs" / f"s{j0['seed']}_zen{format_zen(j0['zen'])}_az{format_zen(j0['az'])}_h{format_height(j0['h'])}_x{_fmt_offset(j0['tel_x'])}_z{_fmt_offset(j0['tel_z'])}.log")
            print(f"DRY-RUN: positional args (slurm $1..$14):")
            print(f"DRY-RUN:   $1  OUTPUT_BASE_DIR  = {output_base}")
            print(f"DRY-RUN:   $2  LOG_FILE         = {log_file_str}")
            print(f"DRY-RUN:   $3  PDG              = {args.pid}")
            print(f"DRY-RUN:   $4  ENERGY           = {args.energy}")
            print(f"DRY-RUN:   $5  ZENITH           = {format_zen(j0['zen'])}")
            print(f"DRY-RUN:   $6  AZIMUTH          = {format_zen(j0['az'])}")
            print(f"DRY-RUN:   $7  INJ_HEIGHT       = {format_height(j0['h'])}")
            print(f"DRY-RUN:   $8  OBS_LEVEL        = {state['config']['obs_level']}")
            print(f"DRY-RUN:   $9  TEL_X            = {j0['tel_x']}")
            print(f"DRY-RUN:   $10 TEL_Y            = {args.tel_y}")
            print(f"DRY-RUN:   $11 TEL_Z            = {j0['tel_z']}")
            print(f"DRY-RUN:   $12 CHERENKOV_RADIUS = {compute_cherenkov_radius(args.energy):.1f}")
            print(f"DRY-RUN:   $13 TEL_RADIUS       = {args.radius}")
            print(f"DRY-RUN:   $14 SEED             = {j0['seed']}")
            print(f"DRY-RUN:   opt HADRON_MODEL    = {state['config']['hadron_model']}")
            print(f"DRY-RUN:   opt RERUN_EVENT     = false")
            print(f"DRY-RUN:   opt --update_time_stamp={str(state['config']['update_time_stamp']).lower()}")
        if bp:
            print(f"DRY-RUN: branch plan: {bp}")
        # Also print the save and next-controller commands
        if jobs:
            work_dir = energy_output_dir(args.pid, args.energy, args.radius) / "adaptive_work" / "cycle_00000"
            config = state["config"]
            n_save = int(config.get("save_batches", 1))
            if n_save > 1:
                slices = split_list(jobs, min(n_save, len(jobs)))
                for idx, batch_jobs in enumerate(slices):
                    if not batch_jobs:
                        continue
                    bchunk = build_chunk_file(work_dir / "batches", batch_jobs, name=f"cycle_chunk_b{idx}.json")
                    bresult = work_dir / "batches" / f"cycle_result_b{idx}.csv"
                    wcmd = build_save_worker_command(bchunk, bresult, config)
                    print(f"DRY-RUN: literal save-batch {idx} sbatch command:")
                    print(f"DRY-RUN: {shlex.join(wcmd)}")
                batch_results = [
                    work_dir / "batches" / f"cycle_result_b{idx}.csv"
                    for idx in range(len(slices))
                ]
                mcmd = build_merge_command(batch_results, config, mode="append")
                print(f"DRY-RUN: literal merge sbatch command:")
                print(f"DRY-RUN: {shlex.join(mcmd)}")
            else:
                chunk_file = build_chunk_file(work_dir, jobs)
                result_csv = work_dir / "cycle_result.csv"
                save_cmd = build_save_command(chunk_file, result_csv, config, mode="append")
                print(f"DRY-RUN: literal save sbatch command:")
                print(f"DRY-RUN: {shlex.join(save_cmd)}")
            next_cmd = build_next_controller_command(args, state_path_for(args.pid, args.energy, args.radius), save_job_id="SAVE_JID_PLACEHOLDER")
            print(f"DRY-RUN: literal next-controller sbatch command:")
            print(f"DRY-RUN: {shlex.join(next_cmd)}")
        return 0

    # ── Normal execution ──
    state_lock = acquire_state_lock(state_file)
    try:
        if state_file.exists():
            state = load_json(state_file)
        else:
            state = build_state_template(args)
            write_json_atomic(state_file, state)
        guard_old_schema(state)

        # Ensure config is current
        state["config"] = default_config(args)
        state["energy_string"] = args.energy
        state["pid"] = args.pid
        state["radius"] = args.radius
        state["tel_y"] = args.tel_y
        state.setdefault("branches", {})
        state.setdefault("shell_history", [])
        state.setdefault("pending", [])
        state.setdefault("status", "active")
        state.setdefault("current_shell", 0)
        state.setdefault("cycle", 0)
        state.setdefault("last_save_job", None)
        ensure_shell_structure(state, int(state["current_shell"]))

        if args.resume:
            print(f"Resume: ignoring {len(state.get('pending', []))} pending entries")
            state["pending"] = []

        # Check if already done
        if state.get("status") == "done":
            print(f"State is 'done' for energy {args.energy}.")
            return 0

        # Ingest CSV data
        if state.get("pending"):
            print(f"Ingesting {len(state['pending'])} pending jobs from previous cycle")
        rows_seen = update_state_from_csvs(state)
        if args.resume:
            print(f"Re-derived from CSVs; rows_seen={rows_seen}")
        normalize_active_branches(state)

        # Advance shell if needed
        maybe_advance_shell(state)
        if state.get("status") == "done":
            print(f"Adaptive scan complete for {args.energy}.")
            write_json_atomic(state_file, state)
            return 0

        # Plan next work
        jobs, branch_plan, _cell_plan = plan_jobs_for_state(state)
        if not jobs:
            state["cycle"] = int(state.get("cycle", 0)) + 1
            write_json_atomic(state_file, state)
            next_cmd = build_next_controller_command(args, state_file, begin_delay=True)
            print(f"Cycle {state['cycle']}: no work; scheduling only next controller with 2min delay.")
            next_jid = submit_with_retry(next_cmd)
            print(f"Next controller job: {next_jid}")
            return 0

        # Check queue
        qproc = subprocess.run(
            ["squeue", "-u", os.environ.get("USER", ""), "-h", "-r"],
            capture_output=True, text=True,
        )
        queued = len([l for l in (qproc.stdout or "").splitlines() if l.strip()])
        n_save = int(state["config"].get("save_batches", 1))
        # Reserve room for the save pipeline + next controller: N worker saves
        # + 1 merge + 1 controller when batching, else save + controller.
        overhead = (n_save + 2) if n_save > 1 else 2
        room = max(0, int(args.max_queued) - queued - overhead)
        if room == 0:
            print(f"Queue at {queued}/{args.max_queued}, no room. Scheduling only next controller with delay.")
            state["cycle"] = int(state.get("cycle", 0)) + 1
            write_json_atomic(state_file, state)
            next_cmd = build_next_controller_command(args, state_file, begin_delay=True)
            next_jid = submit_with_retry(next_cmd)
            print(f"Next controller job: {next_jid}")
            return 0

        planned = jobs
        jobs = jobs[:room]
        if len(jobs) < len(planned):
            print(f"Truncated from {len(planned)} to {len(jobs)} due to max_queued")

        # Build chunk and submit sims
        work_dir = energy_output_dir(args.pid, args.energy, args.radius) / "adaptive_work" / f"cycle_{int(state['cycle']):05d}"
        config = state["config"]
        output_base = energy_output_dir(args.pid, args.energy, args.radius)
        sim_jids = []
        submitted_jobs = []
        for job in jobs:
            cmd = build_sim_command(
                output_base, args.pid, args.energy, args.radius, args.tel_y,
                job["zen"], job["az"], job["h"], job["tel_x"], job["tel_z"], job["seed"],
                config,
            )
            jid = submit_with_retry(cmd)
            sim_jids.append(jid)
            submitted_jobs.append(job)
            print(f"SIM: zen={job['zen']:.1f} h={job['h']} seed={job['seed']} jid={jid}")

        # Track attempts
        apply_attempt_rounds(state, submitted_jobs)

        # Save strategy: batched worker save jobs + one merge job, or the single
        # combined save job when --save-batches <= 1.
        n_save = int(config.get("save_batches", 1))
        if n_save > 1 and jobs:
            slices = split_list(jobs, min(n_save, len(jobs)))
            worker_jids = []
            save_jid = None
            start = 0
            for idx, batch_jobs in enumerate(slices):
                if not batch_jobs:
                    continue
                bchunk = build_chunk_file(work_dir / "batches", batch_jobs, name=f"cycle_chunk_b{idx}.json")
                bresult = work_dir / "batches" / f"cycle_result_b{idx}.csv"
                wcmd = build_save_worker_command(bchunk, bresult, config)
                batch_jids = sim_jids[start:start + len(batch_jobs)]
                start += len(batch_jobs)
                if batch_jids:
                    wrap_idx = wcmd.index("--wrap")
                    wcmd.insert(wrap_idx, f"--dependency=afterany:{':'.join(batch_jids)}")
                wj = submit_with_retry(wcmd)
                worker_jids.append(wj)
                print(f"SAVE-BATCH {idx}: {len(batch_jobs)} events -> {bresult} jid={wj}")
            if worker_jids:
                batch_results = [
                    work_dir / "batches" / f"cycle_result_b{idx}.csv"
                    for idx in range(len(slices))
                ]
                mcmd = build_merge_command(batch_results, config, mode="append")
                wrap_idx = mcmd.index("--wrap")
                mcmd.insert(wrap_idx, f"--dependency=afterany:{':'.join(worker_jids)}")
                save_jid = submit_with_retry(mcmd)
                print(f"SAVE-MERGE: {len(worker_jids)} batch results -> dest, jid={save_jid}")
        else:
            chunk_file = build_chunk_file(work_dir, jobs)
            result_csv = work_dir / "cycle_result.csv"
            save_cmd = build_save_command(chunk_file, result_csv, config, mode="append")
            if sim_jids:
                wrap_idx = save_cmd.index("--wrap")
                save_cmd.insert(wrap_idx, f"--dependency=afterany:{':'.join(sim_jids)}")
            save_jid = submit_with_retry(save_cmd)

        # Submit next controller
        next_cmd = build_next_controller_command(args, state_file, save_job_id=save_jid)
        next_jid = submit_with_retry(next_cmd)

        # Update state
        state["pending"] = [[j["zen"], j["az"], j["h"], j["seed"]] for j in submitted_jobs]
        state["last_save_job"] = save_jid
        state["cycle"] = int(state.get("cycle", 0)) + 1
        write_json_atomic(state_file, state)

        print(f"Cycle {state['cycle']}: {len(submitted_jobs)} sims, save={save_jid}, next={next_jid}")
        return 0
    finally:
        release_state_lock(state_lock)


if __name__ == "__main__":
    raise SystemExit(main())