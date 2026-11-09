# Adaptive Trigger-Region Scanner for Trinity Muon Simulations

This document describes the adaptive scan system that automatically finds the trigger region for Trinity muon air-shower simulations. Instead of running a fixed grid of `(zenith, height)` combinations, the adaptive scanner explores the parameter space intelligently, doubling the number of seeds per cell until a trigger is found or the cell is exhausted, and terminating the scan when a shell produces zero triggers.

---

## 1. Overview

For a given muon energy, the scanner explores the `(zenith, azimuth, height)` parameter space:

- **Zenith / Azimuth** are the **expected arrival** angles of the muon at the detector plane (the input to the magnetic-field correction, see §14). Azimuth is a user-supplied list (`--az-values`); the shipped `adaptive_defaults` carry 7 values, `[270.0, 269, 271, 268, 272, 267, 273]`.
- **Scan order** is `zenith (shells) -> azimuth (list) -> height (N-doubling)`.
- **Telescope position** is fixed at `(tel_x=0, tel_z=0)`; it is relocated later only as a fallback (recentering).
- **Zenith** starts at `zen_start` (default 90.0°) and expands outward in shells, stopping when a shell's branches would exceed `zen_max` (default 92.0°).
- **Height** is scanned from `h_start` upward by `delta_h`, or over an explicit `h_values` list. The shipped default carries a log-spaced `h_values` list, so the arithmetic scan is replaced unless you clear `h_values`.

Every `(zenith, azimuth, height)` cell is converted to a **predicted injection**
by `magnetic_field_muon_correction.back_track_muon` (§14). Cells whose muon
`survival_probability` is below `--min-survival-probability` (shipped default
**`0.5`**; the controller's own argparse default is `0.1`), or whose geometry
cannot be back-tracked, are **skipped** (no jobs submitted). A skipped cell
contributes **no trigger**, so a shell in which *every* cell is skipped produces
zero triggers and ends that energy's scan (§15).

For each cell, the scanner runs simulations with an **N-doubling strategy**: it starts with `n_init` seeds, and if no trigger is found, doubles the number of seeds (`n_init, 2*n_init, 4*n_init, ...`) up to `n_max`. The **shipped defaults are `n_init = n_max = 1`**, i.e. exactly one shower per cell, always seed index `1`; raise `n_init`/`n_max` if you need statistics per cell. The strategy ensures that cells near the trigger boundary get enough statistics to be confident, while cells far from the boundary are quickly exhausted.

The scan is **self-resubmitting**: each controller invocation submits simulation jobs, a save job, and the next controller job via SLURM dependencies. The chain continues until the scan is complete.

---

## 2. File Inventory

| File | Role |
|---|---|
| `Trinity_adaptive_input.py` | Generates the energy list and adaptive defaults as JSON. |
| `adaptive_controller.py` | The state machine controller. One invocation = one cycle. |
| `magnetic_field_muon_correction.py` | `back_track_muon`: arrival angles -> predicted injection + survival probability. |
| `submit_adaptive_scan.sh` | Thin launcher: generates input JSON, accepts any energy string, submits the first controller job. |
| `save_data2csv/save_CARE2csv_chunk_tree.py` | Chunk worker: reads a JSON list of `[zen, az, h, tel_x, tel_z, seed]` tuples, writes one CSV. Supports backward-compatible 5-tuples. |
| `save_data2csv/append_adaptive_rows.py` | Merges cycle-result CSV rows into per-zenith CSVs (append or replace mode). |
| `migrate_tree_remove_injection.sh` | One-off migration: strips the legacy `_ix*_iy*` suffix from old tree dirs. |

---

## 3. State Machine

### 3.1 Definitions

- **Cell** = `(zen, az, h)`. `zen`/`az` are the expected arrival angles; `az` is drawn from the `--az-values` list. `tel_x = tel_z = 0`.
- **Branch** = a fixed `(zen, az)` pair, iterating over `h` (either arithmetic from `h_start` by `delta_h` to `h_max`, or over the explicit `h_values` list).
- **Shell k** (k ≥ 1) = the two zeniths `{zen_start - k*delta_zen, zen_start + k*delta_zen}`, each crossed with the full azimuth list.
  Shell 0 = the single zenith `zen_start` (= 90.0).
- **Trigger in a cell** = ≥1 CSV row for that cell with `file_found==1 AND max_pe >= trigger_pe`. With `min_triggers` > 1, a cell only stops N-doubling when it has **≥ `min_triggers`** such rows.
- **Skipped cell** = the magnetic backtracker found no valid crossing (`skip_reason: invalid_geometry`) or the muon survival is below the threshold (`skip_reason: low_survival`); no jobs are submitted for it.

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
(no-failure) case at exactly `1..target`. A cell is marked `status="failed"`
**only** when it can no longer make forward progress — i.e. `needed` is empty
after the replacement expansion, which happens only if no fresh index is left
to draw from (a seed that never resolves into a CSV row). There is **no longer
a blanket per-target round cap**: cells that are still accumulating successes
toward their doubling `target` (throttled by `--max-queued`) remain `active`
and keep being re-submitted each cycle until they reach `target`.

**Cumulative doubling**: at each target the cell needs at most `target`
distinct indices, drawn from `1..target`. With the shipped defaults
(`n_init=1, n_max=1`) there is exactly one target — `[1..1]` — so each cell runs a
single shower and never doubles. The ladder `[1..n_init]`, `[1..2*n_init]`,
`[1..4*n_init]`, ... up to `n_max` only applies once you raise `n_init`/`n_max`.
All indices are reused across cells. The number of doubling rounds is logged per
cell in the `attempts` dictionary (keyed by target).

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
| `--energy E` | *(required)* | Energy string. Any value is accepted; grid membership is only a NOTE, never an error |
| `--n-init` | 1 | Initial number of seeds per cell |
| `--n-max` | 1 | Maximum number of seeds per cell — with `n_init = n_max = 1` each cell gets exactly **one** shower (seed index `1`) |
| `--min-triggers 1` | 1 | Minimum number of triggered events (`max_pe >= trigger_pe`) required before a cell stops N-doubling. Shell termination still requires ≥1 trigger |
| `--h-start 3000` | 3000 | Starting injection height (m) — ignored when `--h-values` is set |
| `--h-max 100000` | 100000 | Maximum injection height (m) — ignored when `--h-values` is set |
| `--delta-h 10000` | 10000 | Height step (m) — ignored when `--h-values` is set |
| `--h-terminate-above 50000` | 50000 | Height above which a no-trigger cell terminates the branch (only relabels the terminal state when `h_batch >= len(h_values)`, see §5.3) |
| `--h-batch 100` | 100 | Number of consecutive h values processed per cycle. When `>= len(h_values)` the whole height list is a single group (§5.3) |
| `--zen-start 90.0` | 90.0 | Starting zenith (deg) |
| `--zen-max 92.0` | 92.0 | Maximum zenith (deg). The scan stops expanding when a shell's branches would exceed this value. Note there is **no lower bound** (§15) |
| `--h-values 4000,6000,...` | 10 log-spaced heights | Comma-separated list of injection heights (m). When set (non-empty), **replaces** `h_start`/`delta_h`/`h_max` scanning with this explicit list. Shipped list: `[4000, 6000, 8000, 11000, 15000, 21000, 29000, 41000, 57000, 80000]` (deduplicated/sorted/rounded on parse) |
| `--delta-zen 1` | 1 | Zenith step between shells (deg) |
| `--az-values` | 7 values | Arrival azimuths to scan (deg). The launcher passes the `az_values` list from `adaptive_defaults`: `270.0, 269, 271, 268, 272, 267, 273`. `--az` is accepted as an alias |
| `--min-survival-probability 0.5` | 0.5 | Skip a cell when the predicted muon survival probability is below this (see §15). The controller's own argparse default is `0.1` |
| `--trigger-pe 20` | 20.0 | Trigger threshold (PE). Also passed to the save worker as `--Max_PE_cut`, so sub-threshold run dirs are deleted (§7.2) |
| `--max-queued 900` | 900 | Maximum queued jobs |
| `--save-batches` | 100 | Split each cycle's save into N parallel worker jobs + one merge job. `N > 1` means **N save jobs per cycle** — a parallelism knob, *not* "save every N batches" (§5.3) |
| `--tel-x` | 0.0 | Telescope x offset (m) |
| `--tel-z` | 0.0 | Telescope z offset (m) |
| `--resume` | off | Resume from existing state, ignore pending |
| `--rerun-failed` | off | Re-submit failed rows with same seeds |
| `--rerun-detector` | off | Rerun all triggered rows (`max_pe >= trigger_pe`) using only the detector-sim chain (no CORSIKA 8) |
| `--rerun-all` | off | Rerun all triggered rows, including the CORSIKA 8 simulation |
| `--dry-run` | off | Print commands without submitting |
| `--help` | | Show usage |

The launcher also fixes these parameters: `--pid` (= input `pdg`, 13),
`--radius` (= input `tel_radius`, currently **1**), `--tel-y` (= input `tel_y`, 0),
`--obs-level 2944`, and `--hadron-model SIBYLL-2.3d`.

Both `--opt value` and `--opt=value` are accepted.

> **Note on defaults.** The launcher's effective defaults come exclusively from
> `adaptive_defaults` in `Trinity_adaptive_input.py` — **that file is
> authoritative**; the table above mirrors it. As shipped it sets
> `n_init=1, n_max=1, min_triggers=1, h_start=3000, h_max=100000, delta_h=10000,
> h_terminate_above=50000, h_batch=100, zen_start=90.0, zen_max=92.0, delta_zen=1,
> min_survival_probability=0.5, trigger_pe=20.0, max_queued=900, save_batches=100`,
> plus `az_values`, `tel_x`, `tel_z` and a 10-entry `h_values` list. Because
> `h_values` is non-empty the arithmetic scan (`h_start`/`delta_h`/`h_max`) is
> **disabled by default**, and because `h_batch=100` exceeds the list length each
> branch consumes the whole height list in one cycle (see §5.3).

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

### 5.3 Effective defaults and their consequences

The shipped `adaptive_defaults` are not neutral values — several of them change the
shape of the campaign. `Trinity_adaptive_input.py` is authoritative; the table below
records what each value actually does in `adaptive_controller.py`.

| Shipped default | Effective behaviour |
|---|---|
| `n_init = n_max = 1` | Each cell runs exactly **one** shower. `choose_target` returns `None` as soon as `n_success >= 1`, so the cell becomes `exhausted` after a single successful run, and **every cell in every energy uses seed index `1`** (§13). |
| `h_batch = 100` with a 10-entry `h_values` | The whole height list is one group, so a branch is evaluated in a single cycle and is then `complete`. Shell 0 therefore submits 1 zenith × 7 azimuths × 10 heights = **70 sims** in its first cycle; shells 1 and 2 add 2 zeniths × 7 × 10 = **140** each. |
| `h_terminate_above = 50000` | Because all heights sit in the same group, this only relabels the terminal state (`terminated_early` instead of `complete`) — it does not prevent any cell from running. |
| `save_batches = 100` | Each cycle's jobs are split into `min(100, len(jobs))` slices, **each with its own save-worker job**, plus one merge job. For a 70-job cycle that is 70 save jobs (1 event each) + 1 merge. It is a parallelism knob, *not* a "save every N batches" threshold — use ~4-8 for a lighter scheduler footprint. |
| `trigger_pe = 20` | Trigger threshold **and** the save worker's `--Max_PE_cut`: any run with `max_pe < 20` has its whole run directory deleted and its row written with `deleted_low_pe=1, file_found=1` (§7.2). |
| `max_queued = 900` | `room = 900 - queued - (n_save + 2)`, i.e. `900 - queued - 102` (the overhead uses the configured `save_batches`, not the actual slice count). Shell 0 needs 70 sims + 102 overhead ≈ 172 slots, so it does not throttle. |
| `min_survival_probability = 0.5` | The dominant filter of the whole campaign — see §15. |

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
    "n_init": 1, "n_max": 1, "min_triggers": 1,
    "h_start": 3000, "h_max": 100000, "delta_h": 10000,
    "h_terminate_above": 50000, "h_batch": 100,
    "zen_start": 90.0, "zen_max": 92.0, "delta_zen": 1,
    "az_values": [270.0, 269.0, 271.0, 268.0, 272.0, 267.0, 273.0],
    "tel_x": 0.0, "tel_z": 0.0,
    "min_survival_probability": 0.5,
    "trigger_pe": 20.0, "max_queued": 900, "save_batches": 100,
    "obs_level": 2944, "hadron_model": "SIBYLL-2.3d",
    "update_time_stamp": false, "delete_log_on_success": false,
    "h_values": [4000, 6000, 8000, 11000, 15000, 21000, 29000, 41000, 57000, 80000]
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
- **Per-cell seed indices.** Seed index `k` is used at every `(zen, h)` cell; a cell always runs indices `1..target`, where `target` is its current N-doubling target stored in the cell record. There is **no global counter**. The seed reaches CORSIKA **verbatim** (`--seed "${SEED}"` in the slurm script), so cells are not distinguished by the seed at all — independence comes only from their different injection geometry, and with the shipped `n_init = n_max = 1` every cell uses `--seed 1` (see §13). Each cell carries `extra_indices` (default `[]`, populated only by the shortfall path — see §3.2), `seeds` (every seed index ever attempted for the cell, from the CSV), `attempts` (per-target submission round counts, diagnostic only), and `trigger_count`. Both `target` (default `n_init`) and the optional fields default gracefully on read so an older state file still loads.
- **Branch height bookkeeping.** In `h_values` mode a branch advances an index into the list, tracked by `next_h_idx` (not `next_h`). In arithmetic mode it tracks `next_h`, the next arithmetic height.
- **`flock`** the state file for the whole read-modify-write (via a `.lock` sibling).
- **Idempotent**: the controller tolerates being run twice on the same state (e.g., duplicated dependency). Achieved by making the "decide next work" step a pure function of `(state, CSVs)` and clearing `pending` only after successful ingestion.
- **`--resume`**: read the state, ignore `pending`, re-derive work from the CSVs, and restart the chain. Prints what it recovered.
- **No migration**: a state file written with the old global-`seed_counter` schema is rejected at load (`guard_old_schema` raises when `seed_counter` is present) with a clear error and a non-zero exit. Delete the output tree and state file, then relaunch.

---

## 7. CSV Output

### 7.1 Per-(Zenith, Azimuth) CSVs

After each save job, rows are appended to per-cell CSV files, named with both
**expected** (scanned) angles:

```
$BASE_PATH/Muon_pid{pid}_E{E}_R{r}/csv_output/adaptive_pid{pid}_E{E}_R{r}_y{tel_y}_zen{zen:.1f}_az{az:.1f}.csv
```

One cycle can touch two zeniths → rows are split by their expected `scan_zen`
and `scan_az` columns (`zen`/`az` for CSVs predating those columns) and appended
to the matching file. The header is written only when creating the file. `flock`
is used on each cell CSV during append.

> The controller only globs `..._y{tel_y}_zen*_az*.csv`, so per-zenith files from
the old naming scheme are **ignored** rather than double-counted. Those cells
look un-run and are re-simulated (and regenerated with the new columns).

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

### 8.4 `--rerun-detector` / `--rerun-all`

Separate, non-advancing modes (mutually exclusive with each other). They read all
per-zenith CSVs for this energy, select every *triggered* row (`file_found == 1`
**and** `max_pe >= trigger_pe`), and re-submit those exact
`(zen, az, h, tel_x, tel_z, seed)` combos, following with a single save job that
**replaces** the rerun rows in place (mode `replace`). They do **not** touch
`branches`/`current_shell`/`status` and do **not** chain a further controller.

- `--rerun-detector`: passes `RERUN_EVENT=detector_sim` to the slurm script, so the
  existing CORSIKA 8 photon file is **reused** and only the detector-simulation
  chain re-runs (no CORSIKA 8 run).
- `--rerun-all`: passes `RERUN_EVENT=all`, re-running CORSIKA 8 from scratch.

`--dry-run` prints the sim commands, the save command, and the triggered-event count
without submitting anything.

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
| opt | --injection-zenith / --injection-azimuth | predicted injection direction (§14) |
| opt | --injection-x / --injection-y | predicted injection position (m) |
| opt | --injection-z-minus-r | injection z - R_earth (metadata only) |
| opt | --survival-probability | muon survival probability (metadata only) |

The `logs/` directory is created by the controller before sbatch (matching the old driver). When `--rerun-failed` is active, `RERUN_EVENT` is set to `detector_sim`.

The `--injection-*` optionals are produced by the controller from
`magnetic_field_muon_correction.back_track_muon` for each cell (§14). They are
**not** positional args; the positional `$5`/`$6` are the expected arrival angles.

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

> `zen{zen:.1f}`/`az{az:.1f}` are the **expected arrival** angles (scan input);
> `h{h_int}` is the production height. The predicted injection is **not** part of
> the path (it lives in `metadata.yaml`). Legacy `h{h}_ix*_iy*` directories are
> migrated with `migrate_tree_remove_injection.sh <BASE_PATH> --apply`.

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
11. **`bash -n`** on `migrate_tree_remove_injection.sh` and `run_corsika8_trinity_chain_tree.slurm` — no syntax errors.
12. **`py_compile`** on `magnetic_field_muon_correction.py` and `Trinity_adaptive_input.py` — no syntax errors.
13. **Magnetic correction**: `compute_cell_injection` returns a valid injection for downgoing, horizontal and slightly-upgoing zeniths (curved atmosphere); `survival_probability` is `NaN` (cell skipped) only when there is no crossing.
14. **State machine**: branches are keyed by `(zen, az)`; `ensure_shell_structure`, `plan_jobs_for_state`, `current_shell_finished` and `current_shell_any_trigger` verified on a 2-azimuth test shell.
15. **Command generation**: `build_sim_command` emits the `--injection-*` optionals; `build_save_command` no longer passes injection offsets.
16. **Any energy string** accepted by `submit_adaptive_scan.sh` (grid membership is a NOTE, not an error); `Trinity_adaptive_input.py` serialises numpy arrays (`h_values`).

---

## 12. Controlling the Minimum Seeds per Cell (`n_init`)

The minimum number of seeds per cell before moving to the next N-doubling step is controlled by `n_init` in `Trinity_adaptive_input.py`'s `adaptive_defaults` (**currently `1`**). Edit that value in the input file (or pass `--n-init N` to the launcher) to change it. `n_max` (**currently `1`**) similarly bounds the top of the doubling ladder. With the shipped `n_init = n_max = 1` each cell runs exactly one shower and never N-doubles; raise both if you need per-cell statistics. Note that the **launcher** reads these from the input JSON, while the controller's **argparse** defaults (used only when it is invoked directly) are `n_init=1, n_max=64`.

---

## 13. Known Limitations & Notes

- **Shared seed indices across cells**: seed index `k` is reused at every `(zen, h)` cell, so showers at different cells are **not** distinguished by the run-directory `_s{k}` alone. The slurm script passes the seed through **verbatim** (`--seed "${SEED}"` in `run_corsika8_trinity_chain_tree.slurm`) — there is no composite seeding that mixes in zenith or height. Independence across cells therefore comes only from their different injection geometry, and with the shipped `n_init = n_max = 1` every cell in every energy is literally `--seed 1`. Because of this, parallel branches within a shell write into the **same** top-level `_s{k}/` run directories (as sibling `zen*/` subtrees); do not rely on top-level directories being per-`zen`.
- **Monolithic controller**: `adaptive_controller.py` is ~1300 lines. Any future edit should be a targeted patch to a named function, never a full regeneration. If a regeneration ever seems necessary, split the file into modules first.
- **Conda environment**: The controller and save jobs require the `jupyter_env` conda environment (for `uproot`, `pyarrow`, etc.). The `--wrap` commands activate it via `source .../conda.sh && conda activate jupyter_env`.
- **Cherenkov radius**: Currently hard-coded to 15.0 m for all energies. The old driver has a commented-out energy-dependent formula (`max(15.0, 1.5e6 / E)`) that is not active.
- **One energy per chain**: The controller handles exactly one `energy_string`. To scan multiple energies, launch separate chains.
- **`--rerun-failed` is non-advancing**: It does not touch the state machine; it only retries failed rows and replaces them in the CSVs.
- **Default arithmetic scanning disabled**: the shipped `adaptive_defaults` carry a non-empty `h_values` list, so by default the height axis uses that explicit list and `h_start`/`delta_h`/`h_max` are ignored. To restore arithmetic scanning, clear `h_values` in the input JSON (or do not pass `--h-values`).

---

## 14. Magnetic-field muon correction (arrival → injection)

The scan parameters `(zenith, azimuth)` are now the **expected arrival** angles of the
muon at the detector plane, **not** the injection angles. For every `(zenith, azimuth,
height)` cell the controller calls `back_track_muon` (from
`magnetic_field_muon_correction.py`) to predict the **injection** position and direction
that brings the muon to the telescope, and passes that prediction to the CORSIKA 8 job.

### 14.1 Scan order

`zenith (shells) → azimuth (user list) → height (N-doubling)`.

The azimuth list is set with `--az-values` (comma-separated, e.g.
`--az-values 267.0,267.3,270.0`) or the `az_values` key in `Trinity_adaptive_input.py`.
Branches are keyed by `(zen, az)`.

### 14.2 Survival-probability gate

`back_track_muon` also returns a `survival_probability` (muon decay over the
back-tracked path). A cell is **skipped** when

* the geometry cannot be back-tracked (`success == False`, i.e. no crossing within the
  search window) -> `skip_reason: invalid_geometry`, or
* `survival_probability < --min-survival-probability` -> `skip_reason: low_survival`.

The threshold is configurable via `--min-survival-probability` or the
`min_survival_probability` key in `adaptive_defaults`. The **shipped default is
`0.5`** (the controller's own argparse default is `0.1`, used only when it is invoked
directly). Skipped cells are marked `status: skipped` and submit **no jobs**, so they
also contribute no triggers — a shell whose cells are all skipped produces zero
triggers and terminates the scan for that energy. §15 quantifies which
energies/zeniths survive the `0.5` cut.

### 14.3 Tree layout change

The height level no longer carries the injection offsets:

```
<OUTPUT_BASE_DIR>/pdg{PDG}_E{ENERGY}_r{TEL_RADIUS}_s{SEED}/
  zen{ZENITH}/            # EXPECTED arrival zenith (scan input)
    az{AZIMUTH}/          # EXPECTED arrival azimuth (scan input)
      h{INJ_HEIGHT}/      # production height (m)
        x{TEL_X}_y{TEL_Y}_z{TEL_Z}/
          CORSIKA8/ CARE/ GROPT/ CIO/
          metadata.yaml
          correction_report_firstpass.json
```

`migrate_tree_remove_injection.sh <BASE_PATH> [--apply]` strips the legacy
`_ix*_iy*` suffix from existing `h<number>_ix*_iy*` directories (dry-run by default).

### 14.4 CSV columns

`zen`/`az` now hold the **actual arrival** direction of the **centre of mass** of
all particles at the observation level, energy-weighted with
`w_i = kinetic_energy_i * weight_i` (thinning-corrected):
`n_com = sum(w_i n_i) / |sum(w_i n_i)|` from `particles.parquet` `nx/ny/nz`.
`com_concentration` reports `|sum(w_i n_i)| / sum(w_i)` (1.0 = perfectly
collimated) and `com_x_m`/`com_y_m`/`com_z_m` the energy-weighted centroid of the
particle positions. `x_mag`/`y_mag`/`z_mag` still hold the final position of the
highest-kinetic-energy muon (used to centre the telescope), whose direction is
kept in `max_e_muon_zen`/`max_e_muon_az` (the previous meaning of `zen`/`az`).
The expected (scanned) arrival is kept in
`scan_zen`/`scan_az`, and the predicted injection in `injection_zenith_deg`,
`injection_azimuth_deg`, `injection_x`, `injection_y`, `injection_height_m`,
`injection_z_minus_R_m`, plus `survival_probability`.

Consumers that need the scanned cell (completion pre-scans, dedup keys, CARE path
reconstruction for uploads) read `scan_zen`/`scan_az` and only fall back to
`zen`/`az` for legacy CSVs.

The per-cell CSV filename uses both **expected** angles
(`..._zen{ZEN}_az{AZ}.csv`); see §7.1.

### 14.5 CORSIKA 8 injection height

The predicted injection is passed as `--injection-zenith/-azimuth/-x/-y` with
`--injection-height = h_prod` (the rebuilt `muon_cherenkov_3d` interprets
`--injection-height` as the *radial* altitude and derives z). `--injection-z-minus-r`
and `--survival-probability` are carried into `metadata.yaml` for the CSV.

---

## 15. Operating envelope & traps

Measured against the shipped `adaptive_defaults`
(`min_survival_probability = 0.5`, `zen_start = 90.0`, `zen_max = 92.0`,
`h_values` = the 10 log-spaced heights 4-80 km, `obs_level = 2944`).

### 15.1 The survival gate sets the energy floor

`back_track_muon` for pdg `-13`, az 270, obs 2944 m. "kept" = survival ≥ 0.5, i.e.
the cell is actually simulated (heights kept out of the 10 shipped `h_values`):

| E (GeV) | zen 89 | zen 90 | zen 91 | zen 92 |
|---|---|---|---|---|
| 100 | 5/10 (h ≤ 15 km) | 3/10 (h ≤ 8 km) | 0/10 | 0/10 |
| 1e3 | 10/10 | 10/10 | 10/10 | 10/10 |
| 1e4 … 1e8 | 10/10 | 10/10 | 10/10 | 10/10 |

At zen 90 / h 4000 m the transition is sharp:

```
E <= 50 GeV   -> success = False     -> skipped (invalid_geometry)
E = 100 GeV   -> survival = 0.7585   -> KEPT
E = 300 GeV   -> survival = 0.9337   -> KEPT
E = 1 TeV     -> survival = 0.9810   -> KEPT
```

So with the shipped `0.5`:

* Below roughly 60-100 GeV the back-tracker finds no crossing at all (the muon
  decays before reaching the observation level), so **those cells never run**.
* At 100 GeV the zen 91 and zen 92 cells are *entirely* skipped
  (survival 0.42 / 0.20 at h = 4 km).
* From ~1 TeV upward the gate is effectively inactive (20/20 at every zenith).

### 15.2 One fully-skipped shell ends that energy's scan

`maybe_advance_shell` marks the scan `done` when a finished shell produced no
trigger. Because skipped cells submit no jobs, a shell whose cells are *all*
skipped cannot produce a trigger — so a low-energy run can terminate in shell 0
(zen 90) without ever reaching zen 91/92. Of the shipped 112 energies
(0.89 GeV → 8.9e10 GeV), 81 are ≥ 1 TeV; the sub-100 GeV tail contributes almost
nothing.

### 15.3 No lower bound on zenith

`shell_branch_keys` only caps the *upper* side at `zen_max`: `low = zen_start -
shell_index * delta_zen` is added whenever `low <= zen_max`. If triggers keep
occurring, shells continue 90 → 89/91 → 88/92 → 87 → 86 → … → 0 → negative
zeniths. In practice the scan stops when the back-tracker starts failing, but do
not rely on a physical lower limit being enforced.

### 15.4 Cells that never resolve retry indefinitely

With `n_max = 1`, a cell whose save job keeps writing `file_found=0` (e.g. CARE
missing) has `already` filled and `needed` emptied; the controller then appends a
fresh `extra_indices` entry (2, 3, 4, …) one per cycle instead of marking the cell
`failed`. A permanently-broken cell therefore consumes one new sim job per cycle.

### 15.5 Energy grid and naming

`Trinity_adaptive_input.generate_energy_strings()` ships the full hard-coded MCEq
grid (the `1e1…1e3` filter is commented out): **112 values, 0.89 GeV → 8.9e10 GeV**.
`to_1e` uses `int(math.log10(val))`, so sub-GeV energies serialise with an `e0`
exponent and appear in paths as e.g. `E0.891251e0`.

### 15.6 Recommended settings

* Want the low-energy or zen 91/92 end to exist? Lower `min_survival_probability`
  to ~0.05-0.1 (the controller's own default); the shipped `0.5` is what removes them.
* `save_batches = 100` produces ~100 tiny save jobs per cycle — use ~4-8 unless you
  deliberately want that parallelism.
* `n_init = n_max = 1` gives one shower per cell; raise both for per-cell statistics
  (this is also what makes every cell share seed index `1`, see §13).
* Because `h_batch = 100` exceeds `len(h_values)`, `h_terminate_above` cannot
  truncate a branch early; reduce `h_batch` if you want the 50 km cutoff to actually
  stop scanning upward.

