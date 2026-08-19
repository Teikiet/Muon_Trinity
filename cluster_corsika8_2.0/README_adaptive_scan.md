# Adaptive Trigger-Region Scanner for Trinity Muon Simulations

This document describes the adaptive scan system that automatically finds the trigger region for Trinity muon air-shower simulations. Instead of running a fixed grid of `(zenith, height)` combinations, the adaptive scanner explores the parameter space intelligently, doubling the number of seeds per cell until a trigger is found or the cell is exhausted, and terminating the scan when a shell produces zero triggers.

---

## 1. Overview

For a given muon energy, the scanner explores the `(zenith, azimuth, injection_height)` parameter space:

- **Azimuth** is fixed at 270.0°.
- **Telescope position** is fixed at `(tel_x=0, tel_z=0)`.
- **Zenith** starts at 90.0° (vertical) and expands outward in shells, stopping when a shell's branches would exceed `zen_max` (default 92.0°).
- **Height** is scanned from `h_start` (default 4000 m) upward by `delta_h` (default 10000 m). Because the default configuration ships with a non-empty `h_values` list, the **arithmetic scan is replaced by the explicit list** `[3500, 5000, 7000, 9000, 12000, 15000, 20000]` unless you clear `h_values`. The `h_start`/`delta_h`/`h_max` parameters are ignored whenever `h_values` is set.

For each `(zenith, height)` cell, the scanner runs simulations with an **N-doubling strategy**: it starts with `n_init` seeds (default 100), and if no trigger is found, doubles the number of seeds (`n_init, 2·n_init, 4·n_init, ...`) up to `n_max` (default 1000). This ensures that cells near the trigger boundary get enough statistics to be confident, while cells far from the boundary are quickly exhausted.

The scan is **self-resubmitting**: each controller invocation submits simulation jobs, a save job, and the next controller job via SLURM dependencies. The chain continues until the scan is complete.

---

## 2. File Inventory

| File | Role |
|---|---|
| `Trinity_adaptive_input.py` | Generates the energy list and adaptive defaults as JSON. |
| `adaptive_controller.py` | The state machine controller. One invocation = one cycle. |
| `submit_adaptive_scan.sh` | Thin launcher: generates input JSON, validates energy, submits the first controller job. |
| `save_data2csv/save_CARE2csv_chunk_tree.py` | Chunk worker: reads a JSON list of `[zen, az, h, tel_x, tel_z, seed]` tuples, writes one CSV. Supports backward-compatible 5-tuples. |
| `save_data2csv/append_adaptive_rows.py` | Merges cycle-result CSV rows into per-zenith CSVs (append or replace mode). |

---

## 3. State Machine

### 3.1 Definitions

- **Cell** = `(zen, az, h)`. `az` is fixed at 270.0. `tel_x = tel_z = 0`.
- **Branch** = a fixed `zen`, iterating over `h` (either arithmetic from `h_start` by `delta_h` to `h_max`, or over the explicit `h_values` list).
- **Shell k** (k ≥ 1) = the two zeniths `{zen_start - k*delta_zen, zen_start + k*delta_zen}`.
  Shell 0 = the single zenith `zen_start` (= 90.0).
- **Trigger in a cell** = ≥1 CSV row for that cell with `file_found==1 AND max_pe >= trigger_pe`. With `min_triggers` > 1, a cell only stops N-doubling when it has **≥ `min_triggers`** such rows.

### 3.2 Branch Execution (per zenith, independent)

A branch is a sequence of `h` values — an arithmetic series `h_start, h_start+delta_h, ..., <= h_max` (when no `h_values`), or the explicit `h_values` list (when present). Process them in groups of `h_batch` consecutive values (default 100).

For each cell in the current group, run the **N-doubling loop**:

```
targets = [n_init, 2*n_init, 4*n_init, ...] truncated so target <= n_max
for target in targets:
    ensure the cell has `target` SUCCESSFUL runs
    if trigger found in cell: STOP this cell (success)
if no trigger after reaching n_max successful runs: STOP this cell (exhausted)
```

**"Ensure `target` successful runs"**:

**Seed indices are per-cell cumulative.** A cell always draws its seeds from
the contiguous range `1..target`, where `target` is its current N-doubling
target. The seed index `k` is therefore reused at every `(zen, h)` cell. The
downstream slurm script (`run_corsika8_trinity_chain_tree.slurm`) derives
CORSIKA's RNG seed compositely from the seed index together with zenith and
injection height, so reusing an index at different `(zen, h)` still yields
independent showers. The seed index is still passed to the slurm script as
positional argument `$14`, unchanged in meaning.

- `n_success` = rows for this cell with `file_found==1` (in the per-zenith CSV).
- `already` = every seed index already attempted for this cell (all CSV rows
  matching `(zen, az, h)`, **regardless of `file_found`**; persisted in the cell's
  `seeds` field).
- `needed = [k for k in range(1, target+1) if k not in already]` — the untried
  indices in the fixed `1..target` range. Submit **exactly** the indices in
  `needed`.
- After the save step, recount. A failed index is **not** re-derived within the
  fixed range (it stays in `already`/`seeds`); retry it manually with `--rerun-failed`
  (which resubmits the exact index). This is correct behaviour for a transient
  failure but is separate from the shortfall logic below.
- Rows with `file_found==0` remain in the CSV (for `--rerun-failed`) but do
  **not** count toward `target` (they leave `n_success` unchanged).

**Shortfall replacement (`extra_indices`).** Because the index set is fixed at
`1..target`, a permanently-broken index cannot be replaced by a different one
within the range. When the fixed set has been fully attempted but the cell is
still short of `target` successful runs, draw fresh indices from
`target+1, target+2, ...` until `len(needed) + n_success == target`. Record
these in the cell as `extra_indices` so they are **not** re-derived on the next
cycle. This keeps the "N successful runs" contract while keeping the common
(no-failure) case at exactly `1..target`. The 5-round replacement cap is
retained; reaching it while still short marks the cell `status="failed"`.

**Cumulative doubling**: at each target the cell needs at most `target`
distinct indices, drawn from `1..target`. With the current defaults
(`n_init=100, n_max=1000`) and zero failures the per-target ranges are
`[1..100]`, `[1..200]`, `[1..400]`, `[1..800]`, `[1..1000]` (all indices reused
across cells). The number of doubling rounds is logged per cell in the
`attempts` dictionary (keyed by target).

**Branch termination**:
- A cell with no trigger at `n_max` and `h < h_terminate_above` (20000): record it, **continue** to the next `h`.
- A cell with no trigger at `n_max` and `h >= h_terminate_above`: **terminate the branch immediately** (do not scan higher `h`). Mark branch `status="terminated_early"`.
- Reaching the end of the `h` sequence (past `h_max`, or past the end of `h_values`): mark branch `status="complete"`.

### 3.3 Shell Progression

- Shell 0 (`zen = 90.0`) runs first, alone.
- Then shell 1, shell 2, … Within a shell the **two branches run in parallel**: the controller advances both branches' `h`-groups in the same cycle, submitting jobs for both.
- A shell is **joined** only when *both* its branches have `status != "active"`.
- **Shell termination**: if a completed shell produced **zero triggers** across both branches and all `h` cells visited (a truncated branch still counts as a valid verdict), the whole scan is **DONE**. Otherwise proceed to shell k+1.
- Shell 0 producing zero triggers also ends the scan.
- **Zenith limit**: if the next shell's branches would all exceed `zen_max`, the scan is **DONE** (no further expansion).

---

## 4. Controller Cycle

Each invocation of `adaptive_controller.py` does **exactly one** step and exits:

```
1. Load/init  adaptive_state.json  (flock it)
2. If a previous cycle submitted runs (state.pending non-empty):
      -> ingest: read the per-zenith CSVs, update cell counters/triggers,
         clear state.pending, advance the state machine
3. If state.status == "done": print summary, exit 0 (submit nothing)
4. Decide the next work item(s): a list of (zen, az, h, seed) tuples
   (may span 2 branches x h_batch cells x several seeds)
5. If nothing to submit and nothing pending -> state.status="done", exit 0
6. Respect --max-queued: query `squeue -u $USER -h -r | wc -l`.
   Reserve overhead for the save pipeline + next controller: 2 jobs in
   single-save mode (1 save + 1 controller), or (n_save + 2) when batching
   (n_save workers + 1 merge + 1 controller). room = max_queued - queued -
   overhead. If room == 0, submit no sims, and submit only the next controller
   with a short delay (--begin=now+2minutes) so the chain does not die.
   Otherwise TRUNCATE to_submit to room; the remainder is re-derived next cycle
   (state is idempotent, so this is safe).
7. sbatch each sim job  -> collect JIDs
8. Save the cycle's results:
   - If `--save-batches == 1`: sbatch ONE save job
     (`--dependency=afterany:<all sim JIDs>`) that runs the worker then the appender.
   - If `--save-batches N > 1`: split the submitted jobs into N contiguous
     batches; sbatch N **worker-only** save jobs, each depending only on its own
     batch of sim JIDs and writing a temp `cycle_result_b{k}.csv`; then sbatch
     ONE **merge** job (`--dependency=afterany:<all worker JIDs>`) that merges the
     batch result CSVs into the per-zenith destination in a single flocked append.
     The merge JID is the cycle's `save` JID.
9. sbatch the NEXT controller job with --dependency=afterany:<save JID>
   (the merge JID when batching).
10. Record state.pending = to_submit, state.cycle += 1, write state, exit 0
```

**`submit_with_retry`**: retries on `QOSMaxSubmitJobPerUserLimit` (and
`Resource temporarily unavailable` / `Socket timed out`) up to 20 times with
30 s sleep.

---

## 5. CLI Usage

### 5.1 `submit_adaptive_scan.sh`

Thin launcher: generates the input JSON via `Trinity_adaptive_input.py`, validates the requested energy is in the list, then submits the **first** controller job. Everything after that is self-driven.

```bash
bash submit_adaptive_scan.sh --energy E [options]
```

**Options** (all optional except `--energy`). The **Default** column shows the
effective defaults used by the launcher — these come from `adaptive_defaults`
in `Trinity_adaptive_input.py`; the launcher always passes them explicitly, so
they override the controller's own argparse defaults (see §5.2).

| Option | Default (launcher) | Description |
|---|---|---|
| `--energy E` | *(required)* | Energy string, must match an `energy_string` in the input JSON |
| `--n-init 100` | 100 | Initial number of seeds per cell |
| `--n-max 1000` | 1000 | Maximum number of seeds per cell |
| `--min-triggers 1` | 1 | Minimum number of triggered events (`max_pe >= trigger_pe`) required before a cell stops N-doubling. Shell termination still requires ≥1 trigger |
| `--h-start 4000` | 4000 | Starting injection height (m) — ignored when `--h-values` is set |
| `--h-max 100000` | 100000 | Maximum injection height (m) — ignored when `--h-values` is set |
| `--delta-h 10000` | 10000 | Height step (m) — ignored when `--h-values` is set |
| `--h-terminate-above 20000` | 20000 | Height above which a no-trigger cell terminates the branch |
| `--h-batch 100` | 100 | Number of consecutive h values processed per cycle |
| `--zen-start 90.0` | 90.0 | Starting zenith (deg) |
| `--zen-max 92.0` | 92.0 | Maximum zenith (deg). The scan stops expanding when a shell's branches would exceed this value |
| `--h-values 3500,5000,...` | from `adaptive_defaults.h_values` | Comma-separated list of injection heights (m). When set (non-empty), **replaces** `h_start`/`delta_h`/`h_max` scanning with this explicit list. The current default list is `[3500, 5000, 7000, 9000, 12000, 15000, 20000]` (the raw input JSON has a duplicate 15000, deduplicated on parse) |
| `--delta-zen 0.3` | 0.3 | Zenith step between shells (deg) |
| `--az 270.0` | 270.0 | Azimuth (deg) |
| `--trigger-pe 20` | 20.0 | Trigger threshold (PE) |
| `--max-queued 900` | 900 | Maximum queued jobs |
| `--save-batches` | 100 | Split each cycle's save into N parallel worker jobs + one merge job (`--save-batches N` overrides; 1 = single combined save) |
| `--resume` | off | Resume from existing state, ignore pending |
| `--rerun-failed` | off | Re-submit failed rows with same seeds |
| `--dry-run` | off | Print commands without submitting |
| `--help` | | Show usage |

The launcher also fixes these parameters: `--pid` (= input `pdg`, 13),
`--radius` (= input `tel_radius`, currently **1**), `--tel-y` (= input `tel_y`, 0),
`--obs-level 2944`, and `--hadron-model SIBYLL-2.3d`.

Both `--opt value` and `--opt=value` are accepted.

> **Note on defaults.** The launcher's effective defaults come exclusively from
> `adaptive_defaults` in `Trinity_adaptive_input.py`, which currently sets
> `n_init=100, n_max=1000, h_start=4000, h_max=100000, delta_h=10000,
> h_batch=100, save_batches=100` and a non-empty `h_values` list. Because
> `h_values` is non-empty, the arithmetic scan is **disabled by default**.

### 5.2 `adaptive_controller.py` (direct)

The controller can also be invoked directly (e.g., for testing). When invoked
directly **without** the launcher it falls back to its own argparse defaults
(`n_init=1, n_max=64, h_start=3000, h_max=50000, delta_h=200, h_batch=1,
save_batches=1, radius=5`), which differ from the launcher's configured
defaults. Prefer launching through `submit_adaptive_scan.sh` to get the
configured defaults.

```bash
python3 adaptive_controller.py --energy 1.12202e1 --dry-run
```

It accepts all the same options as the launcher plus:
- `--state-file PATH` — override the state file location
- `--pid`, `--radius`, `--tel-y`, `--tel-x`, `--tel-z` — simulation parameters
- `--obs-level 2944`, `--hadron-model SIBYLL-2.3d` — physics parameters
- `--update-time-stamp`, `--delete-log-on-success` — sim job flags

---

## 6. State File

`$BASE_PATH/Muon_pid{pid}_E{E}_R{r}/adaptive_state.json`

```json
{
  "version": 2,
  "energy_string": "1.12202e1",
  "pid": 13,
  "radius": 1,
  "tel_y": 0,
  "config": {
    "n_init": 100, "n_max": 1000, "min_triggers": 1,
    "h_start": 4000, "h_max": 100000, "delta_h": 10000,
    "h_terminate_above": 20000, "h_batch": 100,
    "zen_start": 90.0, "zen_max": 92.0, "delta_zen": 0.3,
    "az": 270.0, "tel_x": 0.0, "tel_z": 0.0,
    "trigger_pe": 20.0, "max_queued": 900, "save_batches": 100,
    "obs_level": 2944, "hadron_model": "SIBYLL-2.3d",
    "update_time_stamp": false, "delete_log_on_success": false,
    "h_values": [3500, 5000, 7000, 9000, 12000, 15000, 20000]
  },
  "status": "active",
  "cycle": 7,
  "current_shell": 2,
  "branches": {
    "89.7": {
      "status": "active",
      "next_h_idx": 2,
      "cells": {
        "7000": {"n_success": 800, "n_failed": 1, "triggered": false, "status": "active", "target": 1000, "attempts": {"800": 1, "1000": 1}, "seeds": [1, 2, 3, 4], "extra_indices": [1001], "trigger_count": 0},
        "9000": {"n_success": 100, "n_failed": 0, "triggered": true, "status": "triggered", "target": 100, "attempts": {"100": 1}, "seeds": [1], "extra_indices": [], "trigger_count": 1}
      },
      "trigger_count": 1,
      "shell": 2
    },
    "90.3": {
      "status": "terminated_early",
      "shell": 2,
      "cells": {}
    }
  },
  "shell_history": [
    {"shell": 0, "zeniths": ["90.0"], "any_trigger": true}
  ],
  "pending": [
    [89.7, 270.0, 7000, 4],
    [89.7, 270.0, 7000, 1001]
  ],
  "last_save_job": "12345678"
}
```

**Requirements:**
- **Per-cell seed indices.** Seed index `k` is used at every `(zen, h)` cell; a cell always runs indices `1..target`, where `target` is its current N-doubling target stored in the cell record. There is **no global counter**. Independence between cells is guaranteed by the slurm script's composite RNG seeding (seed index + zenith + injection height). Each cell carries `extra_indices` (default `[]`, populated only by the shortfall path — see §3.2), `seeds` (every seed index ever attempted for the cell, from the CSV), `attempts` (per-target submission round counts, used for the 5-round cap), and `trigger_count`. Both `target` (default `n_init`) and the optional fields default gracefully on read so an older state file still loads.
- **Branch height bookkeeping.** In `h_values` mode a branch advances an index into the list, tracked by `next_h_idx` (not `next_h`). In arithmetic mode it tracks `next_h`, the next arithmetic height.
- **`flock`** the state file for the whole read-modify-write (via a `.lock` sibling).
- **Idempotent**: the controller tolerates being run twice on the same state (e.g., duplicated dependency). Achieved by making the "decide next work" step a pure function of `(state, CSVs)` and clearing `pending` only after successful ingestion.
- **`--resume`**: read the state, ignore `pending`, re-derive work from the CSVs, and restart the chain. Prints what it recovered.
- **No migration**: a state file written with the old global-`seed_counter` schema is rejected at load (`guard_old_schema` raises when `seed_counter` is present) with a clear error and a non-zero exit. Delete the output tree and state file, then relaunch.

---

## 7. CSV Output

### 7.1 Per-Zenith CSVs

After each save job, rows are appended to per-zenith CSV files:

```
$BASE_PATH/Muon_pid{pid}_E{E}_R{r}/csv_output/adaptive_pid{pid}_E{E}_R{r}_y{tel_y}_zen{zen:.1f}.csv
```

One cycle can touch two zeniths → rows are split by their `zen` column and appended to the matching file. The header is written only when creating the file. `flock` is used on each per-zenith CSV during append.

### 7.2 Row Semantics

| `file_found` | `deleted_low_pe` | `max_pe` | Meaning |
|---|---|---|---|
| 1 | 1 | < cut | Run completed, **no trigger**, run dir deleted |
| 1 | 0 | ≥ cut | Run completed, **TRIGGER**, run dir kept |
| 0 | — | — | Run **FAILED / missing** |

The deletion cut (`--Max_PE_cut`) **must equal** the trigger cut (`--trigger-pe`), so non-triggered runs are deleted every scan and triggered runs are kept.

### 7.3 CSV Fields

The full field list includes: `pid`, `energy_string`, `radius`, `seed`, `zen`, `az`, `height`, `tel_x`, `tel_y`, `tel_z`, `tel_r`, `file_size_MB`, `cph_photon_count`, `cph_max_photons_10ns`, `correction_x_m`, `correction_y_m`, `runtime_seconds`, `max_pe`, `time_at_max_pe_ns`, `avg_pe`, `total_pe`, `n_hit_pixels`, `image_size_pe`, `frac_in_brightest`, `concentration_2`, `pulse_width_ns`, `rise_time_ns`, `time_spread_ns`, `time_gradient`, `peak_to_charge`, `baseline_rms_pe`, `max_muon_energy_GeV`, `r68_m`, `r99_m`, `total_particles`, `muon_component`, `electron_component`, `hadronic_component`, `gamma_component`, `other_component`, `deleted_low_pe`, `file_found`.

---

## 8. Modes

### 8.1 `--dry-run`

Prints every `sbatch` command it *would* run, mutates nothing, submits nothing. Output includes:
- Current shell and active zeniths
- Number of sim jobs planned
- One `CELL` line per cell: `CELL zen=... h=...  target=...  already={...} needed=[...]` listing the per-cell seed indices that would be submitted
- The **literal** sim `sbatch` command for the first job
- All 14 positional arguments with cycle-0 values
- The literal save `sbatch` command (worker + merge commands when `--save-batches > 1`)
- The literal next-controller `sbatch` command

### 8.2 `--resume`

Reads the existing state file, ignores `pending`, re-derives work from the per-zenith CSVs, and restarts the chain. Prints what it recovered.

### 8.3 `--rerun-failed`

Separate, non-advancing mode. When set:
1. Scan all per-zenith CSVs for this energy for rows with `file_found == 0`.
2. Resubmit **those exact `(zen, az, h, seed)`** combos (same seeds — this is a retry, not new statistics).
3. Submit one save job that re-derives their rows and **replaces** them in place in the per-zenith CSVs, keyed on `(zen, az, height, tel_x, tel_z, seed)`.
4. Do **not** touch `branches`, `current_shell`, or `status`. Do **not** chain a further controller.
5. Print a count of rows retried.

---

## 9. Sim Job Arguments (14 positional)

The controller's `build_sim_command` mirrors the call convention in `submit_multiple_continuous_tree.sh` exactly. The slurm script `run_corsika8_trinity_chain_tree.slurm` takes 14 positional args (plus two trailing optional args — the hadron model and the rerun flag):

| # | Slurm var | Cycle-0 value (E=1.12202e1, r=1) |
|---|---|---|
| $1 | OUTPUT_BASE_DIR | `/scratch/general/vast/u1520754/muon_sim_chain_tree_adaptive/Muon_pid13_E1.12202e1_R1` |
| $2 | LOG_FILE | `.../logs/s1_zen90.0_az270.0_h3000_x0_z0.log` |
| $3 | PDG | `13` |
| $4 | ENERGY | `1.12202e1` |
| $5 | ZENITH | `90.0` |
| $6 | AZIMUTH | `270.0` |
| $7 | INJ_HEIGHT | `3000` |
| $8 | OBS_LEVEL | `2944` |
| $9 | TEL_X | `0.0` |
| $10 | TEL_Y | `0` |
| $11 | TEL_Z | `0.0` |
| $12 | CHERENKOV_RADIUS | `15.0` |
| $13 | TEL_RADIUS | `1` |
| $14 | SEED | `1` |
| opt | HADRON_MODEL | `SIBYLL-2.3d` |
| opt | RERUN_EVENT | `false` |
| opt | --update_time_stamp | `false` |

The `logs/` directory is created by the controller before sbatch (matching the old driver). When `--rerun-failed` is active, `RERUN_EVENT` is set to `detector_sim`.

---

## 10. Path Conventions

### 10.1 Simulation Output Tree

```
$BASE_PATH/Muon_pid{pid}_E{energy_str}_R{r}/
  pdg{pid}_E{energy_str}_r{r}_s{seed}/
    zen{zen:.1f}/ az{az:.1f}/ h{h_int}/ x{x}_y{y}_z{z}/
      CARE/cherenkov_hits.root
      CIO/cherenkov_hits.cph
      corsika8_output/particles/particles.parquet
      metadata.yaml
```

`BASE_PATH = /scratch/general/vast/u1520754/muon_sim_chain_tree_adaptive`

> `r` is the telescope radius. The launcher derives `--radius` from the input
> JSON `tel_radius` (currently **1**); the controller's own argparse default is
> 5. `R{r}` in directory names therefore reflects whatever radius the chain was
> launched with.

Because seed indices are reused per cell, one top-level `pdg{pid}_E{energy_str}_r{r}_s{k}/`
directory now holds **multiple `zen*/az*/h*` subtrees** — every cell that used
seed index `k` writes into it. Within a shell the two branches
(`zen = zen_start ± k*delta_zen`) run in parallel and both target the same
top-level `_s{k}/` directories as sibling `zen1/` and `zen2/` subtrees; the
slurm script only creates/populates its own leaf, so concurrent branches do
not overwrite each other.

### 10.2 Controller Work Directory

```
$BASE_PATH/Muon_pid{pid}_E{E}_R{r}/adaptive_work/cycle_{NNNNN}/
  cycle_chunk.json            (single-save mode)
  cycle_result.csv            (single-save mode)
  batches/                    (--save-batches N > 1 mode)
    cycle_chunk_b{k}.json
    cycle_result_b{k}.csv
```

With `--save-batches N > 1`, batch chunk/result files live under `batches/`; the
merge job reads `cycle_result_b{k}.csv` and appends to the per-zenith CSVs.

### 10.3 Log Files

```
$HOME/csv_logs/adaptive_E{E}/cycle_%j.out
$HOME/csv_logs/adaptive_E{E}/cycle_%j.err
$HOME/csv_logs/adaptive_E{E}/save_%j.out
$HOME/csv_logs/adaptive_E{E}/save_%j.err
$HOME/csv_logs/adaptive_E{E}/merge_%j.out   (--save-batches > 1)
$HOME/csv_logs/adaptive_E{E}/merge_%j.err
```

---

## 11. Verification Summary

The following checks were performed and passed against the current codebase:

1. **`python -m py_compile`** on all 4 Python files — no syntax errors.
2. **`bash -n`** on `submit_adaptive_scan.sh` — no syntax errors.
3. **No duplicate `def`s** in `adaptive_controller.py`.
4. **Imports resolve**: all 6 imported names (`_fmt`, `_fmt_angle`, `_fmt_offset`, `make_correction_paths`, `make_path`, `make_run_dir`) exist in `save_CARE2csv_chunk_tree.py`.
5. **6-tuple seed support**: `combo_seed()` present, `--seed default=None`, `combo_seed_value` used everywhere.
6. **14-arg mapping confirmed** against `run_corsika8_trinity_chain_tree.slurm` (plus the two trailing optional args).
7. **Log-file `mkdir`** confirmed before sbatch.
8. **Cherenkov radius** matches old driver (both hard-code 15.0).
9. **`--obs-level` / `--hadron-model`** passed consistently by launcher and chaining function.
10. **State schema `version: 2`** with per-cell `target`/`extra_indices`/`seeds`/`attempts` fields; `guard_old_schema` rejects the removed `seed_counter` schema.

---

## 12. Controlling the Minimum Seeds per Cell (`n_init`)

The minimum number of seeds per cell before moving to the next N-doubling step is controlled by `n_init` in `Trinity_adaptive_input.py`'s `adaptive_defaults` (currently `100`). Edit that value in the input file (or pass `--n-init N` to the launcher) to change it. `n_max` (currently `1000`) similarly bounds the top of the doubling ladder. Note that the **launcher** reads these from the input JSON, while the controller's **argparse** defaults (used only when it is invoked directly) are `n_init=1, n_max=64`.

---

## 13. Known Limitations & Notes

- **Shared seed indices across cells**: seed index `k` is reused at every `(zen, h)` cell, so showers at different cells are **not** distinguished by the run-directory `_s{k}` alone. Independence across cells relies entirely on the slurm script's composite RNG seeding (seed index + zenith + injection height). Because of this, parallel branches within a shell write into the **same** top-level `_s{k}/` run directories (as sibling `zen*/` subtrees); do not rely on top-level directories being per-`zen`.
- **Monolithic controller**: `adaptive_controller.py` is ~1300 lines. Any future edit should be a targeted patch to a named function, never a full regeneration. If a regeneration ever seems necessary, split the file into modules first.
- **Conda environment**: The controller and save jobs require the `jupyter_env` conda environment (for `uproot`, `pyarrow`, etc.). The `--wrap` commands activate it via `source .../conda.sh && conda activate jupyter_env`.
- **Cherenkov radius**: Currently hard-coded to 15.0 m for all energies. The old driver has a commented-out energy-dependent formula (`max(15.0, 1.5e6 / E)`) that is not active.
- **One energy per chain**: The controller handles exactly one `energy_string`. To scan multiple energies, launch separate chains.
- **`--rerun-failed` is non-advancing**: It does not touch the state machine; it only retries failed rows and replaces them in the CSVs.
- **Default arithmetic scanning disabled**: the shipped `adaptive_defaults` carry a non-empty `h_values` list, so by default the height axis uses that explicit list and `h_start`/`delta_h`/`h_max` are ignored. To restore arithmetic scanning, clear `h_values` in the input JSON (or do not pass `--h-values`).
