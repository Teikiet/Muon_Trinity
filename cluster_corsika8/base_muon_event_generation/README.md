# Base Muon Event Generation (Tree Workflow)

This folder contains the Slurm submission and execution scripts used to run the CORSIKA8 -> correction -> Trinity container chain for muon events in a tree-structured output layout.

Primary workflow scripts:
- submit_multiple_continuous_tree.sh
- submit_multiple_continuous_wrapper.slurm
- run_corsika8_trinity_chain_tree.slurm

Related helper:
- ../particle_location_correction.py

## 1. High-level Workflow

The pipeline runs in two layers:

1) Submission layer (many jobs)
- Script: submit_multiple_continuous_tree.sh
- Builds a todo list of (seed, energy, zenith, azimuth, x, z, height) combinations.
- Applies pre-scan skip logic (unless rerun mode is enabled).
- Submits one Slurm job per combination to run_corsika8_trinity_chain_tree.slurm.

2) Per-job execution layer (single combination)
- Script: run_corsika8_trinity_chain_tree.slurm
- Produces/loads Cherenkov eventio data.
- Runs particle location correction (with recentering).
- Optionally reduces base cherenkov_hits_base.dat radius for storage optimization.
- Runs container chain (CIO -> GrOptics -> CARE).
- Flattens output and renames CHASM -> CORSIKA8.

## 2. Output Tree Layout

For each submitted parameter point:

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
- During processing, CHASM exists first and is renamed to CORSIKA8 at the end.
- metadata.yaml is rewritten each run for that node.

## 3. Key Modes in run_corsika8_trinity_chain_tree.slurm

Mode selection is based on TEL_X/TEL_Z and effective rerun state.

Effective rerun is true when either condition is true:
- rerun_event flag is true
- CSV override is triggered: exact combo row has file_found=0 in
  csv_output/scan_care_pid{PDG}_E{ENERGY}_R{TEL_RADIUS}_y{TEL_Y}_s{SEED}.csv

Exact combo match keys used by the per-job script:
- seed, zen, az, height, tel_x, tel_z

A) Base run (TEL_X=0 and TEL_Z=0)
- effective rerun=false: full CORSIKA8 run, then correction, then container chain.
- effective rerun=true: reuses existing base cherenkov_hits_base.dat, skips CORSIKA8.

B) Reuse/offset run (TEL_X!=0 or TEL_Z!=0)
- Reuses base geometry cherenkov_hits_base.dat from x0_y0_z0 node.
- Applies correction for requested offset.

## 4. Submission Script Behavior

Script: submit_multiple_continuous_tree.sh

Main phases:
1) Fast pre-scan
- Reads CSV outputs (file_found column) to mark completed points.
- If --rerun-event is enabled, pre-scan filtering is bypassed and all combos are queued.

2) Per-job effective rerun override
- Even when --rerun-event is not enabled, per-job script can force downstream rerun.
- Trigger: exact CSV row for the submitted combo has file_found=0.
- Effect: job does not early-skip just because CARE output exists.

3) Submission loop
- Computes dynamic CHERENKOV_RADIUS as max(15.0, 1.5e6 / E).
- Enforces queue cap (MAX_QUEUED, POLL_INTERVAL).
- Calls sbatch for each todo point.

Current parameter grid in this script is intentionally minimal:
- zeniths: 89.9
- azimuths: 270.0
- tel_xs: 0
- tel_zs: 0
- heights: 5000
- seeds: 2, 3
- energies: MCEq-derived in [1e3, 2e3] GeV equivalent grid format used by script.

## 5. Flags

## submit_multiple_continuous_tree.sh

Supported flags:
- --rerun-event
- --update_time_stamp
- --reduced_based_run_radius
- --reduced_based_run_radius=<R>

Behavior:
- --rerun-event:
  Forces rerun path in per-job script (skip done checks in submit pre-scan).
  Note: without this flag, per-job script may still force rerun when CSV file_found=0.
- --update_time_stamp:
  Updates mtime on reusable base artifacts when applicable.
- --reduced_based_run_radius[=R]:
  Enables base cherenkov_hits_base.dat reduction pass.
  Default R is 15 meters when no value is given.

## run_corsika8_trinity_chain_tree.slurm optional args

Accepted as optional args after required positional arguments:
- [hadron_model]
- [rerun_event]
- --update_time_stamp[=true|false]
- --reduced_based_run_radius[=R]

## 6. Base DAT Reduction and Safety Logic

When --reduced_based_run_radius is enabled:
- Applies only to base runs (x=0,z=0).
- Validates reduced radius R > 0.
- Validates R <= CHERENKOV_RADIUS for that event.
- Runs correction script on cherenkov_hits_base.dat with telescope-x/y = 0 and radius R.
- Overwrites cherenkov_hits_base.dat with reduced output.

Extra protection for repeat runs:
- If correction_report_firstpass.json already existed at job start,
- and cherenkov_hits_base.dat exists,
- and (rerun_event=true or update_time_stamp=true),
then reduction pass is skipped (assume already reduced).

In that skip path:
- base dat content is reused unchanged.
- if --update_time_stamp is active, both files are touched:
  - cherenkov_hits_base.dat
  - correction_report_firstpass.json

Additional behavior:
- When --reduced_based_run_radius is enabled, script checks whether correction_report_firstpass.json exists.
- If missing at check time, script prints a warning and continues.

## 7. Metadata and Correction Report

metadata.yaml (written per run) includes core run inputs.

Rerun metadata fields:
- rerun_event: effective rerun state used by execution logic
- rerun_event_requested: original input rerun flag value

For base events, after first correction pass, metadata appends:
- recentered_telescope_x_m
- recentered_telescope_y_m

These values are taken from correction_report_firstpass.json generated by particle_location_correction.py with --report-json.
This captures first-pass computed correction center values used for bookkeeping.

correction_report_firstpass.json contains machine-readable correction summary fields, including center_x_m and center_y_m.

## 8. Disk Usage Strategy

To limit storage growth:
- Persistent large duplicate files (.orig, .preserved) are not kept.
- Temporary backups are created only during active correction and deleted afterward.
- Final cleanup removes known backup artifacts if they remain.
- Optional base-radius reduction shrinks cherenkov_hits_base.dat footprint.

## 9. Common Run Commands

From your home path:

Submit standard rerun + timestamp update + default reduced base radius (15 m):

bash $HOME/Muon_Trinity/cluster_corsika8/base_muon_event_generation/submit_multiple_continuous_tree.sh --rerun-event --update_time_stamp --reduced_based_run_radius

Submit with custom reduced base radius (example 12 m):

bash $HOME/Muon_Trinity/cluster_corsika8/base_muon_event_generation/submit_multiple_continuous_tree.sh --rerun-event --update_time_stamp --reduced_based_run_radius=12

Wrapper job (already configured):

sbatch $HOME/Muon_Trinity/cluster_corsika8/base_muon_event_generation/submit_multiple_continuous_wrapper.slurm

## 10. Failure Cases to Know

- Base dat missing in reuse/rerun mode:
  Job exits with error and asks for completed base run first.

- Invalid reduced radius:
  Job exits if reduced radius is not numeric, <= 0, or > CHERENKOV_RADIUS.

- Missing correction_report_firstpass.json with reduced-radius enabled:
  Job logs a warning and continues.

- Correction failure:
  Job restores working dat from temporary backup and exits non-zero.

- Container failure:
  Job exits with container return code after logging final status.

## 11. Quick File Map

- submit_multiple_continuous_tree.sh
  Batch discovery + queue-aware job submission.

- submit_multiple_continuous_wrapper.slurm
  Wrapper Slurm job that launches submit_multiple_continuous_tree.sh.

- run_corsika8_trinity_chain_tree.slurm
  Per-parameter full processing chain.

- ../particle_location_correction.py
  Eventio hotspot finding, radius filtering, optional recentering, JSON reporting.
