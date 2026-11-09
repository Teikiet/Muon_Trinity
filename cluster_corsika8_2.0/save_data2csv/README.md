# save_data2csv Workflow

This folder contains batch scripts and Python workers that convert Trinity CARE outputs into CSV tables.

Main workflows:
- Tree CSV export (current): chunked SLURM arrays over zen/az/height
- Legacy CSV export (older flat-path layout)
- Triggered-event resampling CSV export

## 1. Script Map

Tree workflow (recommended):
- submit_save_care2csv_tree.sh
- save_CARE2csv_chunk_tree.py
- merge_csv_chunks.py

Legacy workflow:
- submit_save_care2csv.sh
- save_CARE2csv_chunk.py
- merge_csv_chunks.py

Triggered-event resampling:
- submit_trigger_resample.sh (chunked array approach)
- build_trigger_resample_chunk.py (worker used by chunked approach)
- merge_trigger_resample_csv.py
- run_build_trigger_csv.slurm (single-job loop approach)
- build_trigger_resample_csv.py (non-chunked builder)

## 2. Data Layout Assumptions

Tree workflow expects simulation output under (adaptive scan base path):

/scratch/general/vast/u1520754/muon_sim_chain_tree_adaptive/
  Muon_pid{PDG}_E{ENERGY_STR}_R{RADIUS}/
    pdg{PDG}_E{ENERGY_STR}_r{RADIUS}_s{SEED}/
      zen{ZEN}/az{AZ}/h{H}/x{X}_y{Y}_z{Z}/
        CORSIKA8/cherenkov_hits.dat
        CORSIKA8/particles/particles.parquet
        CARE/cherenkov_hits.root
        metadata.yaml
        correction_report_firstpass.json

Here `zen{ZEN}`/`az{AZ}` are the EXPECTED arrival angles (scan input) and `h{H}`
is the production height. The legacy `h{H}_ix*_iy*` level is no longer used; run
`migrate_tree_remove_injection.sh <BASE_PATH> --apply` on old trees.

CSV output goes to:

.../Muon_pid{PDG}_E{ENERGY_STR}_R{RADIUS}/csv_output/

### 2.1 CSV schema (adaptive per-cell files)

Per-(zenith, azimuth) CSVs are named
`adaptive_pid{PDG}_E{ENERGY}_R{RADIUS}_y{TEL_Y}_zen{ZEN}_az{AZ}.csv`, where
`ZEN` and `AZ` are the **expected** (scanned) arrival angles at one decimal,
formatted with `_fmt_angle` (`scan_zen`/`scan_az`, falling back to `zen`/`az`
for CSVs written before those columns existed). Columns of interest:

| Column | Meaning |
|---|---|
| `scan_zen`, `scan_az` | EXPECTED arrival angles (scan input; used for the filename and the dedup key) |
| `zen`, `az` | ACTUAL arrival direction of the **centre of mass** of all particles at the observation level, energy-weighted: `w_i = kinetic_energy_i * weight_i` (thinning-corrected), `n_com = sum(w_i n_i) / |sum(w_i n_i)|` (from `CORSIKA8/particles/particles.parquet`, `nx/ny/nz`) |
| `com_concentration` | `|sum(w_i n_i)| / sum(w_i)` for that average; 1.0 = perfectly collimated |
| `max_e_muon_zen`, `max_e_muon_az` | Legacy `zen`/`az`: ACTUAL arrival direction of the single highest-kinetic-energy muon |
| `com_x_m`, `com_y_m`, `com_z_m` | Energy-weighted centroid of all particle positions (telescope-centred `x/y/z`) |
| `x_mag`, `y_mag`, `z_mag` | ACTUAL arrival position of the highest-kinetic-energy muon (telescope-centred `x/y/z`), used to centre the telescope |
| `injection_zenith_deg`, `injection_azimuth_deg` | PREDICTED injection direction (from `back_track_muon`) |
| `injection_x`, `injection_y` | PREDICTED injection position (global x/y, m) |
| `injection_height_m` | production height (radial altitude, m) |
| `injection_z_minus_R_m` | injection z - R_earth (m) |
| `survival_probability` | muon survival probability from `back_track_muon` |
| `correction_x_m`, `correction_y_m` | telescope recentering offset (from the correction report) |

The injection columns are read from each run's `metadata.yaml`; the arrival
columns from `particles.parquet`.

> **Note (COM direction change).** `zen`/`az` used to hold the direction of the
> highest-kinetic-energy muon. They now hold the energy-weighted centre-of-mass
> direction of *all* particles at the observation level. Consumers that need the
> *scanned* cell (path reconstruction, completion pre-scans, dedup keys) must use
> `scan_zen`/`scan_az` and fall back to `zen`/`az` only for legacy CSVs — this is
> implemented in `merge_csv_chunks.py`, `upload_triggered_care_to_drive.py`,
> `submit_save_care2csv_tree.sh` and
> `../base_muon_event_generation/submit_multiple_continuous_tree.sh`.
> The legacy direction is still available as `max_e_muon_zen`/`max_e_muon_az`.
> Pre-existing merged CSVs are mixed-version; regenerate them with the worker
> (no re-simulation needed) to get consistent values.

## 3. Tree CSV Export (Primary)

Entry point:
- submit_save_care2csv_tree.sh

What it does:
1. Builds MCEq energy strings in [1e3, 1e6].
2. Creates chunk JSON files of parameter combinations.
3. Submits SLURM array jobs in batches with queue throttling and retry.
4. Runs merge job after chunk arrays finish.

Default grid in submit_save_care2csv_tree.sh:
- zeniths: 87.0 to 89.7 step 0.3, plus 89.9
- azimuths: 267.0 to 273.0 step 0.3
- tel_x: 0
- tel_z: 0
- heights: 5000 to 50000 (100 points)
- seeds: 1

Key controls:
- CHUNK_SIZE=100
- MAX_JOBS=500
- MAX_SUBMIT_RETRIES=20
- SLEEP_SEC=30

Output file name:
- scan_care_pid{PDG}_E{ENERGY_STR}_R{RADIUS}_y{TEL_Y}_s{SEED}.csv

### file_found meaning in tree worker

In save_CARE2csv_chunk_tree.py, file_found is 1 only when both are true for that row:
- CARE metrics were successfully read from CARE/cherenkov_hits.root
- correction_report_firstpass.json exists and center_x_m/center_y_m were read

Otherwise file_found is 0.

## 4. Legacy CSV Export (Older Layout)

Entry point:
- submit_save_care2csv.sh

Differences from tree workflow:
- Uses base path muon_sim_chain (not muon_sim_chain_tree)
- Uses energy-magnitude format E1e{ENERGY_MAG}
- Worker is save_CARE2csv_chunk.py
- file_found is based on CARE read success only (no correction report requirement)

Use this only if your data are in the old Tilt_* path structure.

## 5. Triggered-Event Resampling CSV

Goal:
- Start from base merged CSV
- Keep events with max_pe >= threshold
- For each triggered base event, probe many (tel_x, tel_z) resample positions
- Build *_trigger.csv

There are two implementations.

A) Chunked array approach (recommended for scale)
- submit_trigger_resample.sh
- Worker: build_trigger_resample_chunk.py
- Merge: merge_trigger_resample_csv.py

Behavior:
1. Reads original merged CSV:
   scan_care_pid{PDG}_E{ENERGY_STR}_R{RADIUS}_y{TEL_Y}_s{SEED}.csv
2. Filters triggered rows by PE_THRESHOLD (default 20).
3. Splits triggered rows into chunk JSON files.
4. Runs array workers to evaluate tel_x/tel_z grids.
5. Merges into:
   scan_care_pid{PDG}_E{ENERGY_STR}_R{RADIUS}_y{TEL_Y}_s{SEED}_trigger.csv

B) Single-job loop approach
- run_build_trigger_csv.slurm calls build_trigger_resample_csv.py in nested loops.

Note:
- submit_trigger_resample.sh uses energy range [1e3, 4e3]
- run_build_trigger_csv.slurm uses energy range [1e3, 1e4]
Adjust as needed for consistency.

## 6. Environment Requirements

Expected environment:
- Conda env: jupyter_env
- Python packages: uproot, numpy, MCEq, crflux
- SLURM access to owner-guest / kingspeak-guest

Scripts source:
- /uufs/chpc.utah.edu/common/home/u1520754/miniconda3/etc/profile.d/conda.sh

Logs:
- $HOME/csv_logs

## 7. Common Commands

Tree export:

bash $HOME/Muon_Trinity/cluster_corsika8_2.0/save_data2csv/submit_save_care2csv_tree.sh

Legacy export:

bash $HOME/Muon_Trinity/cluster_corsika8_2.0/save_data2csv/submit_save_care2csv.sh

Triggered resample (chunked arrays):

bash $HOME/Muon_Trinity/cluster_corsika8_2.0/save_data2csv/submit_trigger_resample.sh

Triggered resample (single SLURM job):

sbatch $HOME/Muon_Trinity/cluster_corsika8_2.0/save_data2csv/run_build_trigger_csv.slurm

## 8. Failure Modes and Notes

- Missing CARE file/tree/branches:
  Worker writes zeros and file_found=0 for affected rows.

- Missing correction report (tree worker):
  correction_x_m and correction_y_m are blank and file_found=0.

- Low-PE deletion (adaptive path):
  the adaptive controller invokes this worker with `--Max_PE_cut <trigger_pe>`, so
  any run whose `max_pe` falls below the trigger threshold has its **whole run
  directory deleted** and its row written with `deleted_low_pe=1, file_found=1`.
  See `README_adaptive_scan.md` §7.2 for the `file_found`/`deleted_low_pe`/`max_pe`
  truth table, and §5.3 for why the cut must equal `--trigger-pe`.

- Partial array submission failure:
  submit_save_care2csv_tree.sh skips merge when not all chunk arrays submit.

- Merge behavior:
  merge scripts concatenate result_*.csv and keep only the first header.

- Re-runs:
  chunk directories are removed and rebuilt for each seed/energy in submit scripts.
