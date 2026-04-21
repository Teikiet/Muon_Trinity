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

Tree workflow expects simulation output under:

/scratch/general/vast/u1520754/muon_sim_chain_tree/
  Muon_pid{PDG}_E{ENERGY_STR}_R{RADIUS}/
    pdg{PDG}_E{ENERGY_STR}_r{RADIUS}_s{SEED}/
      zen{ZEN}/az{AZ}/h{H}/x{X}_y{Y}_z{Z}/
        CARE/cherenkov_hits.root
        correction_report_firstpass.json

CSV output goes to:

.../Muon_pid{PDG}_E{ENERGY_STR}_R{RADIUS}/csv_output/

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

bash $HOME/Muon_Trinity/cluster_corsika8/save_data2csv/submit_save_care2csv_tree.sh

Legacy export:

bash $HOME/Muon_Trinity/cluster_corsika8/save_data2csv/submit_save_care2csv.sh

Triggered resample (chunked arrays):

bash $HOME/Muon_Trinity/cluster_corsika8/save_data2csv/submit_trigger_resample.sh

Triggered resample (single SLURM job):

sbatch $HOME/Muon_Trinity/cluster_corsika8/save_data2csv/run_build_trigger_csv.slurm

## 8. Failure Modes and Notes

- Missing CARE file/tree/branches:
  Worker writes zeros and file_found=0 for affected rows.

- Missing correction report (tree worker):
  correction_x_m and correction_y_m are blank and file_found=0.

- Partial array submission failure:
  submit_save_care2csv_tree.sh skips merge when not all chunk arrays submit.

- Merge behavior:
  merge scripts concatenate result_*.csv and keep only the first header.

- Re-runs:
  chunk directories are removed and rebuilt for each seed/energy in submit scripts.
