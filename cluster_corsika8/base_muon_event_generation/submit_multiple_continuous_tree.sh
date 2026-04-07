#!/bin/bash
# =============================================================================
# submit_multiple_continuous_tree.sh
#
# Usage:
#   bash submit_multiple_continuous_tree.sh
# =============================================================================

# =============================================================================
# FIXED PARAMETERS
# =============================================================================
compute_cherenkov_radius() {
    local ENERGY="$1"
    python3 -c "
import math
E = float('${ENERGY}')
r = max(15.0, 1.5e6 / E)
print(f'{r:.1f}')
"
}
PDG=13
TEL_Y=0
#CHERENKOV_RADIUS=15
TEL_RADIUS=5
OBS_LEVEL=2944
HADRON_MODEL="SIBYLL-2.3d"

SLURM_SCRIPT="$HOME/Muon_Trinity/cluster_corsika8/base_muon_event_generation/run_corsika8_trinity_chain_tree.slurm"

BASE_DIR="/scratch/general/vast/u1520754/muon_sim_chain_tree"

# =============================================================================
# QUEUE CONTROL
# =============================================================================
MAX_QUEUED=500
POLL_INTERVAL=60

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
SEEDS=(1)

# =============================================================================
# PHASE 1: FAST PRE-SCAN (all energies at once)
# =============================================================================

TODO_FILE=$(mktemp /tmp/todo_combos_XXXXXX.txt)

echo "============================================================"
echo "Phase 1: Fast pre-scan using CSV output files..."
echo "============================================================"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate jupyter_env
python3 - "${BASE_DIR}" "${TODO_FILE}" "${SEEDS[*]}" \
          "${PDG}" "${TEL_Y}" "${TEL_RADIUS}" <<'PYEOF'
import sys, os, csv
from itertools import product
import numpy as np
from MCEq.core import MCEqRun
import MCEq.config as config
import math
config.kernel_config= 'MKL'
config.e_min = 1000 #0.160 # GeV
import crflux.models as pm
config.integrator= 'euler'

base_dir    = sys.argv[1]
todo_file   = sys.argv[2]
seeds       = sys.argv[3].split()
pdg         = sys.argv[4]
tel_y       = sys.argv[5]
tel_radius  = sys.argv[6]

# --- All parameter grids ---
def to_1e(val):
    exp = int(math.log10(val))
    coeff = val / 10**exp
    return f"{coeff:g}e{exp}"

# Initialize MCEq with custom atmosphere and Frisco Peak location
mceq = MCEqRun(
    # interaction interaction model
    interaction_model='SIBYLL23C',
    # Primary cosmic ray model
    primary_model=(pm.GlobalSplineFitBeta, None),
    # Set to 0° for horizontal muons
    theta_deg=0.,
    # Use custom atmosphere and geometry
    density_model=("CORSIKA", ('USStd', None)),
    #density_model=("MSIS00_IC", ('FriscoPeak', 'January')),
)
E = mceq.e_grid
E_max = 2e3
E_min = 1e3
E = E[E <= E_max]
E = E[E >= E_min]
energies = [to_1e(e) for e in E]
zeniths  = [f"{z:.1f}" for z in np.arange(87.0, 90, 0.3)]
zeniths.append("89.9")
azimuths = [f"{a:.1f}" for a in np.arange(267.0, 273.3, 0.3)]
tel_xs   = ["0"]
tel_zs   = ["0"]
heights  = [str(int(h)) for h in np.linspace(5000, 50000, 100)]

def key_canon(seed, energy, zen, az, height, tel_x, tel_z):
    return (
        str(int(float(seed))),
        energy,
        f"{float(zen):.1f}",
        f"{float(az):.1f}",
        str(int(round(float(height)))),
        str(int(round(float(tel_x)))) if abs(float(tel_x) - round(float(tel_x))) < 1e-9 else f"{float(tel_x):g}",
        str(int(round(float(tel_z)))) if abs(float(tel_z) - round(float(tel_z))) < 1e-9 else f"{float(tel_z):g}",
    )

completed = set()

# --- Check CSVs for each energy ---
for energy in energies: 
    csv_dir = os.path.join(base_dir, f"Muon_pid{pdg}_E{energy}_R{tel_radius}", "csv_output")

    for seed in seeds:
        csv_path = os.path.join(
            csv_dir,
            f"scan_care_pid{pdg}_E{energy}_R{tel_radius}_y{tel_y}_s{seed}.csv"
        )
        if not os.path.exists(csv_path):
            continue

        with open(csv_path, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row.get("file_found", "0") == "1":
                    key = key_canon(
                        row["seed"], energy, row["zen"], row["az"],
                        row["height"], row["tel_x"], row["tel_z"],
                    )
                    completed.add(key)

        print(f"  Loaded {csv_path}: {len(completed)} total completed", file=sys.stderr)

# --- Check tree structure for each energy ---
tree_completed = 0
for energy in energies:
    output_base = os.path.join(base_dir, f"Muon_pid{pdg}_E{energy}_R{tel_radius}")
    for seed, zen, az, tx, tz, h in product(seeds, zeniths, azimuths, tel_xs, tel_zs, heights):
        key = key_canon(seed, energy, zen, az, h, tx, tz)
        if key in completed:
            continue
        tree_path = os.path.join(
            output_base,
            f"pdg{pdg}_E{energy}_r{tel_radius}_s{seed}",
            f"zen{zen}", f"az{az}", f"h{h}",
            f"x{tx}_y{tel_y}_z{tz}",
            "CARE", "cherenkov_hits.root"
        )
        if os.path.exists(tree_path):
            completed.add(key)
            tree_completed += 1

if tree_completed > 0:
    print(f"  Additional completed from tree structure: {tree_completed}", file=sys.stderr)
print(f"  Total completed combos: {len(completed)}", file=sys.stderr)

# --- Write todo list (now includes energy column) ---
total = 0
done  = 0
todo  = 0

with open(todo_file, "w") as f:
    for energy, seed, zen, az, tx, tz, h in product(energies, seeds, zeniths, azimuths, tel_xs, tel_zs, heights):
        total += 1
        key = key_canon(seed, energy, zen, az, h, tx, tz)
        if key in completed:
            done += 1
        else:
            f.write(f"{seed}\t{energy}\t{zen}\t{az}\t{tx}\t{tz}\t{h}\n")
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
# PHASE 2: SUBMISSION LOOP (reads energy from todo file)
# =============================================================================

SUBMITTED=0
FAILED=0

echo ""
echo "============================================================"
echo "Phase 2: Submitting ${TODO_COUNT} jobs (tree structure)"
echo "  PDG=${PDG}  CHERENKOV_RADIUS=dynamic  TEL_RADIUS=${TEL_RADIUS}"
echo "  Max queued jobs : ${MAX_QUEUED}"
echo "============================================================"
echo ""

while IFS=$'\t' read -r SEED ENERGY ZENITH AZIMUTH TEL_X TEL_Z INJ_HEIGHT; do
    # Compute energy-dependent Cherenkov radius
    CHERENKOV_RADIUS=$(compute_cherenkov_radius "${ENERGY}")
    # Per-energy output directory
    OUTPUT_BASE_DIR="${BASE_DIR}/Muon_pid${PDG}_E${ENERGY}_R${TEL_RADIUS}"
    LOG_DIR="${OUTPUT_BASE_DIR}/logs"
    mkdir -p "${LOG_DIR}"

    LOG_FILE="${LOG_DIR}/s${SEED}_zen${ZENITH}_az${AZIMUTH}_h${INJ_HEIGHT}_x${TEL_X}_z${TEL_Z}.log"

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
        SUBMITTED=$((SUBMITTED + 1))

        if (( SUBMITTED % 50 == 0 )); then
            NJOBS=$(squeue -u "$USER" -h 2>/dev/null | wc -l)
            echo "  [$(date +%H:%M:%S)] Submitted: ${SUBMITTED}/${TODO_COUNT} | Queue: ${NJOBS} | E=${ENERGY} s=${SEED} zen=${ZENITH} az=${AZIMUTH} h=${INJ_HEIGHT}"
        fi
    else
        echo "    WARNING: sbatch failed: ${SBATCH_OUTPUT}"
        FAILED=$((FAILED + 1))

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

rm -f "${TODO_FILE}"

echo ""
echo "============================================================"
echo "Submission complete."
echo "  Todo (from pre-scan) : ${TODO_COUNT}"
echo "  Successfully submitted: ${SUBMITTED}"
echo "  Failed                : ${FAILED}"
echo "============================================================"
