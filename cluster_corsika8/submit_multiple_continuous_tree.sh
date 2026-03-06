#!/bin/bash
# =============================================================================
# submit_multiple_continuous_tree.sh
#
# Phase 1: Fast pre-scan using Python — checks CSV output files for completed runs.
# Phase 2: Continuously submit from the todo list, throttled to MAX_QUEUED.
#
# Tree structure:
#   <OUTPUT_BASE_DIR>/pdg${PDG}_E${ENERGY}_r${TEL_RADIUS}_s${SEED}/
#     zen${ZENITH}/az${AZIMUTH}/h${INJ_HEIGHT}/x${TEL_X}_y${TEL_Y}_z${TEL_Z}/
#
# Usage:
#   bash submit_multiple_continuous_tree.sh
# =============================================================================

# =============================================================================
# FIXED PARAMETERS (same for all jobs)
# =============================================================================
PDG=13                   # muon-
ENERGY=1e5              # primary energy in eV
TEL_Y=0                  # telescope y position (fixed)
CHERENKOV_RADIUS=15     # large sampling area for CORSIKA 8 photon collection (m)
TEL_RADIUS=5             # actual telescope aperture radius (m)
OBS_LEVEL=2944           # observation level in meters
HADRON_MODEL="SIBYLL-2.3d"

# Energy magnitude for CSV path (extract from ENERGY)
ENERGY_MAG=5             # corresponds to 1e5

SLURM_SCRIPT="$HOME/Muon_Trinity/cluster_corsika8/run_corsika8_trinity_chain_tree.slurm"

# Parent output directory (all jobs share this)
BASE_DIR="/scratch/general/vast/u1520754/muon_sim_chain_tree"
OUTPUT_DIR_NAME="Muon_pid${PDG}_E${ENERGY}_R${TEL_RADIUS}"
OUTPUT_BASE_DIR="${BASE_DIR}/${OUTPUT_DIR_NAME}"
# Log directory (one log per job)
LOG_DIR="${OUTPUT_BASE_DIR}/logs"
mkdir -p "${LOG_DIR}"

# CSV output directory (from the CARE2csv pipeline)
CSV_DIR="${BASE_DIR}/Muon_pid${PDG}_E1e${ENERGY_MAG}_R${TEL_RADIUS}/csv_output"

# =============================================================================
# QUEUE CONTROL
# =============================================================================
MAX_QUEUED=500    # maximum jobs in queue at any time
POLL_INTERVAL=60  # seconds between queue checks when full

wait_for_queue_room() {
    while true; do
        NJOBS=$(squeue -u "$USER" -h 2>/dev/null | wc -l)
        if [ "${NJOBS}" -lt "${MAX_QUEUED}" ]; then
            return
        fi
        echo "    [$(date +%H:%M:%S)] Queue at ${NJOBS}/${MAX_QUEUED}. Waiting ${POLL_INTERVAL}s..."
        sleep "${POLL_INTERVAL}"
    done
}

# =============================================================================
# VARIABLE PARAMETERS
# =============================================================================

SEEDS=(1) #($(seq 1 26))

# =============================================================================
# PHASE 1: FAST PRE-SCAN — check CSV files for completed runs
# =============================================================================

TODO_FILE=$(mktemp /tmp/todo_combos_XXXXXX.txt)

echo "============================================================"
echo "Phase 1: Fast pre-scan using CSV output files..."
echo "  CSV directory: ${CSV_DIR}"
echo "============================================================"

python3 - "${OUTPUT_BASE_DIR}" "${TODO_FILE}" "${SEEDS[*]}" \
          "${PDG}" "${ENERGY}" "${TEL_Y}" "${TEL_RADIUS}" \
          "${CSV_DIR}" "${ENERGY_MAG}" <<'PYEOF'
import sys, os, csv
from itertools import product
import numpy as np

output_base = sys.argv[1]
todo_file   = sys.argv[2]
seeds       = sys.argv[3].split()
pdg         = sys.argv[4]
energy      = sys.argv[5]
tel_y       = sys.argv[6]
tel_radius  = sys.argv[7]
csv_dir     = sys.argv[8]
energy_mag  = sys.argv[9]

zeniths  = [f"{z:.1f}" for z in np.arange(87.0, 90, 0.3)]  # 87.0, 87.3, ..., 89.9
zeniths.append("89.9")  # Ensure 89.9 is included
azimuths = [f"{a:.1f}" for a in np.arange(267.0, 273.1, 0.3)]  # 268.0, 268.3, ..., 272.0
tel_xs   = ["0"] #"-5", "-2", "-1", "-0.7", "-0.5", "-0.2", "0", "0.2", "0.5", "0.7", "1", "2", "5"
tel_zs   = ["0"]
heights  = [str(int(h)) for h in np.linspace(3000, 100000, 100)]


def key_canon(seed, zen, az, height, tel_x, tel_z):
    return (
        str(int(float(seed))),
        f"{float(zen):.1f}",
        f"{float(az):.1f}",
        str(int(round(float(height)))),
        str(int(round(float(tel_x)))) if abs(float(tel_x) - round(float(tel_x))) < 1e-9 else f"{float(tel_x):g}",
        str(int(round(float(tel_z)))) if abs(float(tel_z) - round(float(tel_z))) < 1e-9 else f"{float(tel_z):g}",
    )

# ─── Build set of completed combos from CSV files ───
completed = set()

for seed in seeds:
    csv_path = os.path.join(
        csv_dir,
        f"scan_care_pid{pdg}_E1e{energy_mag}_R{tel_radius}_y{tel_y}_s{seed}.csv"
    )
    if not os.path.exists(csv_path):
        print(f"  WARNING: CSV not found for seed {seed}: {csv_path}", file=sys.stderr)
        continue

    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("file_found", "0") == "1":
                key = key_canon(
                    row["seed"], row["zen"], row["az"],
                    row["height"], row["tel_x"], row["tel_z"],
                )
                completed.add(key)

    print(f"  Loaded {csv_path}: found entries contributing to {len(completed)} total completed", file=sys.stderr)

print(f"  Total completed combos from CSVs: {len(completed)}", file=sys.stderr)

# ─── Also check tree structure for completed runs (CARE output exists) ───
tree_completed = 0
for seed, zen, az, tx, tz, h in product(seeds, zeniths, azimuths, tel_xs, tel_zs, heights):
    key = key_canon(seed, zen, az, h, tx, tz)
    if key in completed:
        continue
    # Check tree path for CARE output
    tree_path = os.path.join(
        output_base,
        f"pdg{pdg}_E{energy}_r{tel_radius}_s{seed}",
        f"zen{zen}",
        f"az{az}",
        f"h{h}",
        f"x{tx}_y{tel_y}_z{tz}",
        "CARE",
        "cherenkov_hits.root"
    )
    if os.path.exists(tree_path):
        completed.add(key)
        tree_completed += 1

if tree_completed > 0:
    print(f"  Additional completed from tree structure: {tree_completed}", file=sys.stderr)
print(f"  Total completed combos: {len(completed)}", file=sys.stderr)

# ─── Check each combo against completed set ───
total = 0
done  = 0
todo  = 0

with open(todo_file, "w") as f:
    for seed, zen, az, tx, tz, h in product(seeds, zeniths, azimuths, tel_xs, tel_zs, heights):
        total += 1

        key = key_canon(seed, zen, az, h, tx, tz)

        if key in completed:
            done += 1
        else:
            f.write(f"{seed}\t{zen}\t{az}\t{tx}\t{tz}\t{h}\n")
            todo += 1

        if total % 10000 == 0:
            print(f"  Checked {total} | done {done} | todo {todo}", file=sys.stderr)

print(f"  Total combos : {total}")
print(f"  Already done : {done}")
print(f"  To submit    : {todo}")
PYEOF

TODO_COUNT=$(wc -l < "${TODO_FILE}")

if [ "${TODO_COUNT}" -eq 0 ]; then
    echo ""
    echo "Nothing to submit — all runs already exist!"
    rm -f "${TODO_FILE}"
    exit 0
fi

# =============================================================================
# PHASE 2: SUBMISSION LOOP — read todo list, continuously fill queue
# =============================================================================

SUBMITTED=0
FAILED=0

echo ""
echo "============================================================"
echo "Phase 2: Submitting ${TODO_COUNT} jobs (tree structure)"
echo "  Output base dir : ${OUTPUT_BASE_DIR}"
echo "  SLURM script    : ${SLURM_SCRIPT}"
echo "  PDG=${PDG}  ENERGY=${ENERGY}  CHERENKOV_RADIUS=${CHERENKOV_RADIUS}  TEL_RADIUS=${TEL_RADIUS}"
echo "  Max queued jobs : ${MAX_QUEUED}"
echo "  Tree layout     : pdg_E_r_s/zen/az/h/x_y_z/"
echo "============================================================"
echo ""

while IFS=$'\t' read -r SEED ZENITH AZIMUTH TEL_X TEL_Z INJ_HEIGHT; do
    # Log file uses tree-like naming for easy identification
    LOG_FILE="${LOG_DIR}/s${SEED}_zen${ZENITH}_az${AZIMUTH}_h${INJ_HEIGHT}_x${TEL_X}_z${TEL_Z}.log"

    # Wait until there's room in the queue
    wait_for_queue_room

    # Submit the job
    SBATCH_OUTPUT=$(sbatch "${SLURM_SCRIPT}" \
        "${OUTPUT_BASE_DIR}" \
        "${LOG_FILE}" \
        "${PDG}" \
        "${ENERGY}" \
        "${ZENITH}" \
        "${AZIMUTH}" \
        "${INJ_HEIGHT}" \
        "${OBS_LEVEL}" \
        "${TEL_X}" \
        "${TEL_Y}" \
        "${TEL_Z}" \
        "${CHERENKOV_RADIUS}" \
        "${TEL_RADIUS}" \
        "${SEED}" \
        "${HADRON_MODEL}" 2>&1)

    if [ $? -eq 0 ]; then
        JOBID=$(echo "${SBATCH_OUTPUT}" | awk '{print $NF}')
        SUBMITTED=$((SUBMITTED + 1))

        # Print progress every 50 jobs
        if (( SUBMITTED % 50 == 0 )); then
            NJOBS=$(squeue -u "$USER" -h 2>/dev/null | wc -l)
            echo "  [$(date +%H:%M:%S)] Submitted: ${SUBMITTED}/${TODO_COUNT} | Queue: ${NJOBS} | s=${SEED} zen=${ZENITH} az=${AZIMUTH} h=${INJ_HEIGHT}"
        fi
    else
        echo "    WARNING: sbatch failed: ${SBATCH_OUTPUT}"
        FAILED=$((FAILED + 1))

        # Retry after waiting if it was a queue limit error
        if echo "${SBATCH_OUTPUT}" | grep -q "QOSMaxSubmitJobPerUserLimit"; then
            echo "    Queue limit hit. Waiting 120s and retrying..."
            sleep 120
            wait_for_queue_room

            SBATCH_OUTPUT=$(sbatch "${SLURM_SCRIPT}" \
                "${OUTPUT_BASE_DIR}" \
                "${LOG_FILE}" \
                "${PDG}" \
                "${ENERGY}" \
                "${ZENITH}" \
                "${AZIMUTH}" \
                "${INJ_HEIGHT}" \
                "${OBS_LEVEL}" \
                "${TEL_X}" \
                "${TEL_Y}" \
                "${TEL_Z}" \
                "${CHERENKOV_RADIUS}" \
                "${TEL_RADIUS}" \
                "${SEED}" \
                "${HADRON_MODEL}" 2>&1)

            if [ $? -eq 0 ]; then
                JOBID=$(echo "${SBATCH_OUTPUT}" | awk '{print $NF}')
                echo "    -> Retry succeeded. Job ID: ${JOBID}"
                SUBMITTED=$((SUBMITTED + 1))
                FAILED=$((FAILED - 1))
            else
                echo "    WARNING: Retry also failed: ${SBATCH_OUTPUT}"
            fi
        fi
    fi

done < "${TODO_FILE}"

# Cleanup
rm -f "${TODO_FILE}"

# =============================================================================
# SUMMARY
# =============================================================================

echo ""
echo "============================================================"
echo "Submission complete."
echo "  Todo (from pre-scan) : ${TODO_COUNT}"
echo "  Successfully submitted: ${SUBMITTED}"
echo "  Failed                : ${FAILED}"
echo ""
echo "All outputs will be under (tree structure):"
echo "  ${OUTPUT_BASE_DIR}/pdg${PDG}_E${ENERGY}_r${TEL_RADIUS}_s<SEED>/zen<Z>/az<A>/h<H>/x<X>_y${TEL_Y}_z<Z>/"
echo ""
echo "Per-job SLURM logs (stdout/stderr) are in:"
echo "  /scratch/general/vast/u1520754/log/corsika8_trinity_<jobid>.out/err"
echo "Per-job chain logs are in:"
echo "  ${LOG_DIR}/"
echo "============================================================"
