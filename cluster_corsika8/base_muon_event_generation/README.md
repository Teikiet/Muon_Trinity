# Base Muon Event Generation (Tree Workflow)

This folder contains Slurm scripts for the tree-structured muon workflow:

1. CORSIKA8 eventio generation (or base-dat reuse)
2. particle location correction
3. Trinity container chain (CIO -> GrOptics -> CARE)
4. final flatten/rename to the persistent tree layout

Primary workflow scripts:
- submit_multiple_continuous_tree.sh
- submit_multiple_continuous_wrapper.slurm
- run_corsika8_trinity_chain_tree.slurm

Related helper:
- ../particle_location_correction.py

## 1. High-level Workflow

The pipeline runs in two layers.

1) Submission layer (many jobs)
- Script: submit_multiple_continuous_tree.sh
- Builds a todo list over (seed, energy, zenith, azimuth, x, z, height).
- Applies CSV-based completion filtering unless --rerun-event is enabled.
- Submits one Slurm job per todo combo to run_corsika8_trinity_chain_tree.slurm.

2) Per-job execution layer (single combo)
- Script: run_corsika8_trinity_chain_tree.slurm
- Resolves run mode and effective rerun state.
- Performs cleanup of previous outputs.
- Reuses base dat or runs full CORSIKA8 as needed.
- Applies correction and optional base-dat radius reduction.
- Runs container chain and final output flattening/rename.

## 2. Output Tree Layout

For each submitted point:

<OUTPUT_BASE_DIR>/pdg{PDG}_E{ENERGY}_r{TEL_RADIUS}_s{SEED}/
  zen{ZENITH}/
    az{AZIMUTH}/
      h{INJ_HEIGHT}/
        x{TEL_X}_y{TEL_Y}_z{TEL_Z}/
          CORSIKA8/
          CARE/
          GROPT/
          CIO/
          metadata.yaml
          correction_report_firstpass.json

Notes:
- During execution, CHASM/ is used and renamed to CORSIKA8/ at the end.
- metadata.yaml is rewritten each run for the node.

## 3. Per-job Run Logic (run_corsika8_trinity_chain_tree.slurm)

### 3.1 Effective rerun

EFFECTIVE_RERUN_EVENT starts from rerun_event input, then may be forced true by CSV override.

CSV override trigger:
- CSV row match on seed, zen, az, height, tel_x, tel_z
- row file_found == 0 in:
  csv_output/scan_care_pid{PDG}_E{ENERGY}_R{TEL_RADIUS}_y{TEL_Y}_s{SEED}.csv

### 3.2 Early skip checks

Before cleanup:
- Base run (x=0,z=0) and effective rerun=false:
  skip entire job if both CARE/cherenkov_hits.root and base CORSIKA8/cherenkov_hits_base.dat exist.
- Non-base run and effective rerun=false:
  skip if CARE/cherenkov_hits.root exists.

No early skip is applied when effective rerun=true.

### 3.3 Cleanup behavior

If SIM_DIR already exists, script removes these directories when present:
- CHASM
- CARE
- GROPT
- CIO
- corsika8_output

CORSIKA8 removal rule:
- Remove CORSIKA8 unless (effective rerun=true AND base run=true).

Also removes known artifacts:
- corsika8_tables.dat
- cherenkov_hits_base.dat.preserved
- CORSIKA8/cherenkov_hits.dat.orig
- CHASM/cherenkov_hits.dat.orig

### 3.4 Mode behavior in SECTION 3

A) Non-base run (TEL_X!=0 or TEL_Z!=0)
- Always reuse base geometry cherenkov_hits_base.dat from x0_y0_z0 candidates.
- If missing, job exits with error.

B) Base run (TEL_X=0 and TEL_Z=0)
- effective rerun=true:
  try to reuse base cherenkov_hits_base.dat from candidate locations.
  if found: copy to CHASM working paths.
  if missing: print warning and fallback to FULL RUN (set EFFECTIVE_RERUN_EVENT=false).
- effective rerun=false:
  run full CORSIKA8 and produce new base dat.

Important:
- The startup "Run mode" banner is printed before SECTION 3 fallback.
- metadata.yaml rerun_event field is written before fallback; it reflects pre-fallback effective rerun state.

### 3.5 Correction and reduction

- Runs particle_location_correction.py on CHASM cherenkov_hits.dat with recentering.
- Writes correction_report_firstpass.json.
- For base runs, appends recentered_telescope_x_m and recentered_telescope_y_m to metadata.yaml.

If --reduced_based_run_radius is enabled:
- Applies only on base runs.
- Requires reduced radius > 0 and <= CHERENKOV_RADIUS.
- Operates on CHASM/cherenkov_hits_base.dat.
- Skip reduction when all are true:
  - correction_report_firstpass.json existed before SECTION 4
  - base dat exists
  - (effective rerun=true or update_time_stamp=true)
- In skip path with update_time_stamp=true, touches base dat and correction report.

### 3.6 Container and finalization

- Temporarily moves cherenkov_hits_base.dat out of CHASM before container run, restores afterward.
- Runs apptainer chain script.
- Flattens nested Tilt_* under CHASM/CARE/GROPT/CIO.
- Renames CHASM -> CORSIKA8 (removing existing CORSIKA8 first if present).
- Removes temporary backup artifacts.

## 4. Submission Script Logic (submit_multiple_continuous_tree.sh)

### 4.1 Flags supported

- --rerun-event
- --update_time_stamp (also accepts --update-time-stamp)
- --reduced_based_run_radius
- --reduced_based_run_radius=<R>

### 4.2 Pre-scan phase

- Builds completed set only from CSV rows with file_found=1.
- If --rerun-event is true, skips completion filtering and submits all combos.
- Uses CSV only (no filesystem fallback check).

### 4.3 Parameter grid currently in script

- seeds: 2, 3
- energies: MCEq e_grid values filtered to 1e3 <= E <= 5e3, formatted as coeffeexp (example: 4.46684e3)
- zeniths: 87.0 to 89.7 in 0.3 steps, plus 89.9
- azimuths: 267.0 to 273.0 in 0.3 steps
- tel_xs: 0
- tel_zs: 0
- heights: 100 values linearly spaced from 5000 to 50000 (integer string)

### 4.4 Submission phase

- Computes CHERENKOV_RADIUS per combo as max(15.0, 1.5e6 / E).
- Enforces queue ceiling with MAX_QUEUED and POLL_INTERVAL.
- Submits one sbatch per todo line.
- On QOSMaxSubmitJobPerUserLimit, waits and retries once for that job.

## 5. Optional Args for run_corsika8_trinity_chain_tree.slurm

After required positional args, accepted optionals are:
- [hadron_model]
- [rerun_event]
- --update_time_stamp[=true|false]
- --hadron_model=<model>
- --rerun_event=<true|false>
- --reduced_based_run_radius[=R]

## 6. Metadata and Reports

metadata.yaml contains run input/state fields including:
- rerun_event
- rerun_event_requested

For base events, after correction, metadata appends:
- recentered_telescope_x_m
- recentered_telescope_y_m

correction_report_firstpass.json is produced by particle_location_correction.py with --report-json and contains correction summary fields including center_x_m and center_y_m.

## 7. Disk Usage Strategy

- Avoids keeping persistent .orig/.preserved duplicates.
- Uses temporary backup files during correction and removes them.
- Final cleanup removes residual backup artifacts.
- Optional base-radius reduction can reduce base dat size.

## 8. Common Commands

Submit rerun + timestamp update + default reduced base radius (15 m):

bash $HOME/Muon_Trinity/cluster_corsika8/base_muon_event_generation/submit_multiple_continuous_tree.sh --rerun-event --update_time_stamp --reduced_based_run_radius

Submit with custom reduced base radius (example 12 m):

bash $HOME/Muon_Trinity/cluster_corsika8/base_muon_event_generation/submit_multiple_continuous_tree.sh --rerun-event --update_time_stamp --reduced_based_run_radius=12

Submit wrapper job (currently configured to pass --update_time_stamp --reduced_based_run_radius --rerun-event):

sbatch $HOME/Muon_Trinity/cluster_corsika8/base_muon_event_generation/submit_multiple_continuous_wrapper.slurm

## 9. Failure/Warning Cases

- Base dat missing in non-base reuse mode:
  job exits (base run must exist first).

- Base dat missing in base rerun mode:
  job does not exit; it warns and falls back to full CORSIKA8 for that base event.

- Invalid reduced radius:
  job exits if non-numeric, <= 0, or > CHERENKOV_RADIUS.

- Missing correction_report_firstpass.json with reduced radius enabled:
  warning only; script continues.

- Correction failure:
  restores working dat from temporary backup and exits non-zero.

- Container failure:
  exits with container return code after final logging.

## 10. Quick File Map

- submit_multiple_continuous_tree.sh
  batch discovery + queue-aware submission.

- submit_multiple_continuous_wrapper.slurm
  wrapper Slurm job that launches submit_multiple_continuous_tree.sh.

- run_corsika8_trinity_chain_tree.slurm
  per-combo processing and tree finalization.

- ../particle_location_correction.py
  eventio correction/filtering/recentering + JSON reporting.
