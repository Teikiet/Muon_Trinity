#!/bin/bash
# =============================================================================
# submit_multiple_continuous_tree.sh
#
# Usage:
#   bash submit_multiple_continuous_tree.sh [--rerun-event[=detector_sim|all]] [--update_time_stamp] [--reduced_based_run_radius[=R]] [--only-energy E] [--only-seed S]
# =============================================================================

# =============================================================================
# FIXED PARAMETERS 13:-muon, -13:+muon, 22:photon, 2212:proton
# =============================================================================
compute_cherenkov_radius() {
    local ENERGY="$1"
    python3 -c "
import math
E = float('${ENERGY}')
r = 15.0 #max(15.0, 1.5e6 / E)
print(f'{r:.1f}')
"
}
PDG=""
TEL_Y=""
#CHERENKOV_RADIUS=15
TEL_RADIUS=""
OBS_LEVEL=2944
HADRON_MODEL="SIBYLL-2.3d"
RERUN_EVENT="false"
UPDATE_TIME_STAMP="false"
REDUCED_BASED_RUN_RADIUS_ENABLED="false"
REDUCED_BASED_RUN_RADIUS="15"
ONLY_ENERGY=""
ONLY_SEED=""

# Optional flags
while [[ $# -gt 0 ]]; do
    case "$1" in
        --rerun-event|--rerun_event)
            if [[ $# -gt 1 && ! "$2" =~ ^-- ]]; then
                RERUN_EVENT="$2"
                shift 2
            else
                RERUN_EVENT="detector_sim"
                shift
            fi
            ;;
        --rerun-event=*|--rerun_event=*)
            RERUN_EVENT="${1#*=}"
            shift
            ;;
        --update_time_stamp|--update-time-stamp)
            UPDATE_TIME_STAMP="true"
            shift
            ;;
        --reduced_based_run_radius)
            REDUCED_BASED_RUN_RADIUS_ENABLED="true"
            REDUCED_BASED_RUN_RADIUS="15"
            if [[ $# -gt 1 && ! "$2" =~ ^-- ]]; then
                REDUCED_BASED_RUN_RADIUS="$2"
                shift 2
            else
                shift
            fi
            ;;
        --reduced_based_run_radius=*)
            REDUCED_BASED_RUN_RADIUS_ENABLED="true"
            REDUCED_BASED_RUN_RADIUS="${1#*=}"
            shift
            ;;
        --only-energy|--energy)
            if [[ $# -gt 1 && ! "$2" =~ ^-- ]]; then
                ONLY_ENERGY="$2"
                shift 2
            else
                echo "ERROR: --only-energy requires a value"
                exit 1
            fi
            ;;
        --only-energy=*|--energy=*)
            ONLY_ENERGY="${1#*=}"
            shift
            ;;
        --only-seed|--seed)
            if [[ $# -gt 1 && ! "$2" =~ ^-- ]]; then
                ONLY_SEED="$2"
                shift 2
            else
                echo "ERROR: --only-seed requires a value"
                exit 1
            fi
            ;;
        --only-seed=*|--seed=*)
            ONLY_SEED="${1#*=}"
            shift
            ;;
        --help|-h)
            echo "Usage: bash submit_multiple_continuous_tree.sh [--rerun-event[=detector_sim|all]] [--update_time_stamp] [--reduced_based_run_radius[=R]] [--only-energy E] [--only-seed S]"
            echo "  --rerun-event[=detector_sim|all]   detector_sim: reuse CORSIKA8 output; all: rerun CORSIKA8 + downstream"
            echo "  --update_time_stamp   Touch existing base cherenkov_hits_base.dat to refresh mtime"
            echo "  --reduced_based_run_radius[=R]   Reduce base cherenkov_hits_base.dat with radius R (default: 15 m)"
            echo "  --only-energy E   Run only a single energy string"
            echo "  --only-seed S     Run only a single seed"
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option '$1'"
            echo "Usage: bash submit_multiple_continuous_tree.sh [--rerun-event[=detector_sim|all]] [--update_time_stamp] [--reduced_based_run_radius[=R]] [--only-energy E] [--only-seed S]"
            exit 1
            ;;
    esac
done

case "${RERUN_EVENT,,}" in
    true|1|yes|y)
        RERUN_EVENT="detector_sim"
        ;;
esac

if [[ ! "${REDUCED_BASED_RUN_RADIUS}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "ERROR: --reduced_based_run_radius must be a positive number (meters). Got '${REDUCED_BASED_RUN_RADIUS}'"
    exit 1
fi

if awk "BEGIN {exit !(${REDUCED_BASED_RUN_RADIUS} > 0)}"; then
    :
else
    echo "ERROR: --reduced_based_run_radius must be > 0. Got '${REDUCED_BASED_RUN_RADIUS}'"
    exit 1
fi

SLURM_SCRIPT="$HOME/Muon_Trinity/cluster_corsika8/base_muon_event_generation/run_corsika8_trinity_chain_tree.slurm"
SIM_INPUT_SCRIPT="$HOME/Muon_Trinity/cluster_corsika8/Trinity_sim_imput.py"

BASE_DIR="/scratch/general/vast/u1520754/muon_sim_chain_tree"

# =============================================================================
# QUEUE CONTROL
# =============================================================================
MAX_QUEUED=900
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
SEEDS=()

# =============================================================================
# PHASE 1: FAST PRE-SCAN (all energies at once)
# =============================================================================

TODO_FILE=$(mktemp /tmp/todo_combos_XXXXXX.txt)

echo "============================================================"
echo "Phase 1: Fast pre-scan using CSV output files..."
echo "============================================================"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate jupyter_env
SIM_INPUT_JSON=$(mktemp /tmp/trinity_sim_input_XXXXXX.json)
python3 "${SIM_INPUT_SCRIPT}" --output "${SIM_INPUT_JSON}"

read -r PDG TEL_Y TEL_RADIUS SEEDS_CSV <<< "$(python3 - "${SIM_INPUT_JSON}" <<'PYEOF'
import json
import sys

with open(sys.argv[1], "r", encoding="utf-8") as f:
    data = json.load(f)

pdg_vals = data.get("pdg", [])
pdg = pdg_vals[0] if isinstance(pdg_vals, list) and pdg_vals else data.get("pdg", "")
tel_y = data.get("tel_y", "")
tel_radius = data.get("tel_radius", "")
seeds = data.get("seeds", [])
seed_csv = ",".join(str(s) for s in seeds)

print(f"{pdg} {tel_y} {tel_radius} {seed_csv}")
PYEOF
)"

IFS="," read -ra SEEDS <<< "${SEEDS_CSV}"

ENERGY_FILE=$(mktemp)
python3 - "${SIM_INPUT_JSON}" <<'PYEOF' > "${ENERGY_FILE}"
import json
import sys

with open(sys.argv[1], "r", encoding="utf-8") as f:
    data = json.load(f)

for energy in data.get("energy_strings", []):
    print(energy)
PYEOF

mapfile -t ENERGY_STRS < "${ENERGY_FILE}"
rm -f "${ENERGY_FILE}"
if [ -n "${ONLY_ENERGY}" ]; then
    FOUND_ENERGY=false
    for ENERGY_STR in "${ENERGY_STRS[@]}"; do
        if [ "${ENERGY_STR}" = "${ONLY_ENERGY}" ]; then
            FOUND_ENERGY=true
            break
        fi
    done
    if [ "${FOUND_ENERGY}" = false ]; then
        echo "ERROR: --only-energy '${ONLY_ENERGY}' not found in Trinity_sim_imput.py energy list."
        exit 1
    fi
fi
if [ -n "${ONLY_SEED}" ]; then
    FOUND_SEED=false
    for SEED in "${SEEDS[@]}"; do
        if [ "${SEED}" = "${ONLY_SEED}" ]; then
            FOUND_SEED=true
            break
        fi
    done
    if [ "${FOUND_SEED}" = false ]; then
        echo "ERROR: --only-seed '${ONLY_SEED}' not found in Trinity_sim_imput.py seed list."
        exit 1
    fi
fi
FILTER_ENERGY="${ONLY_ENERGY}" FILTER_SEED="${ONLY_SEED}" SIM_INPUT_JSON="${SIM_INPUT_JSON}" python3 - "${BASE_DIR}" "${TODO_FILE}" "${SEEDS[*]}" \
          "${PDG}" "${TEL_Y}" "${TEL_RADIUS}" "${RERUN_EVENT}" <<'PYEOF'
import sys, os, csv
import json
from itertools import product
import numpy as np

base_dir    = sys.argv[1]
todo_file   = sys.argv[2]
seeds       = sys.argv[3].split()
pdg         = sys.argv[4]
tel_y       = sys.argv[5]
tel_radius  = sys.argv[6]
rerun_event = sys.argv[7].strip().lower() in (
    "true", "1", "yes", "y", "detector_sim", "detector-sim", "all"
)

sim_input_path = os.environ.get("SIM_INPUT_JSON")
if not sim_input_path:
    raise SystemExit("SIM_INPUT_JSON is not set")

with open(sim_input_path, "r", encoding="utf-8") as f:
    sim = json.load(f)

energies = sim.get("energy_strings", [])
geom = sim.get("geometry", {})
zeniths = geom.get("zeniths_deg", [])
azimuths = geom.get("azimuths_deg", [])
tel_xs = geom.get("tel_xs_m", [])
tel_zs = geom.get("tel_zs_m", [])
heights = geom.get("heights_m", [])
filter_energy = os.environ.get("FILTER_ENERGY", "").strip()
filter_seed = os.environ.get("FILTER_SEED", "").strip()

def normalize_seed(val):
    try:
        return str(int(float(val)))
    except Exception:
        return str(val)

def key_canon(seed, energy, zen, az, height, tel_x, tel_z):
    return (
        normalize_seed(seed),
        energy,
        f"{float(zen):.1f}",
        f"{float(az):.1f}",
        str(int(round(float(height)))),
        str(int(round(float(tel_x)))) if abs(float(tel_x) - round(float(tel_x))) < 1e-9 else f"{float(tel_x):g}",
        str(int(round(float(tel_z)))) if abs(float(tel_z) - round(float(tel_z))) < 1e-9 else f"{float(tel_z):g}",
    )

if filter_energy:
    energies = [e for e in energies if e == filter_energy]
if filter_seed:
    filter_seed_norm = normalize_seed(filter_seed)
    seeds = [s for s in seeds if normalize_seed(s) == filter_seed_norm]

if filter_energy or filter_seed:
    print(
        f"  Filters: energy={filter_energy or 'all'} seed={filter_seed or 'all'}",
        file=sys.stderr,
    )

if not energies:
    raise SystemExit("No energies left after filtering")
if not seeds:
    raise SystemExit("No seeds left after filtering")

completed = set()

if rerun_event:
    print("  RERUN mode: skipping completed-output pre-scan filters", file=sys.stderr)
else:
    # --- Check CSVs for each energy ---
    for energy in energies:
        csv_dir = os.path.join(base_dir, f"Muon_pid{pdg}_E{energy}_R{tel_radius}", "csv_output")
        energy_csv_rows = 0
        energy_missing_file_found = 0

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
                    energy_csv_rows += 1
                    file_found_val = row.get("file_found", "0").strip()
                    if file_found_val == "1":
                        key = key_canon(
                            row["seed"], energy, row["zen"], row["az"],
                            row["height"], row["tel_x"], row["tel_z"],
                        )
                        completed.add(key)
                    else:
                        energy_missing_file_found += 1

            print(f"  Loaded {csv_path}: {len(completed)} total completed", file=sys.stderr)

        print(
            f"  Energy {energy}: CSV rows={energy_csv_rows}, file_found=0 rows={energy_missing_file_found}",
            file=sys.stderr,
        )

    print("  Tree fallback disabled: trusting CSV file_found only", file=sys.stderr)
    print(f"  Total completed combos (CSV only): {len(completed)}", file=sys.stderr)

# --- Write todo list (now includes energy column) ---
total = 0
done  = 0
todo  = 0

with open(todo_file, "w") as f:
    for energy, seed, zen, az, tx, tz, h in product(energies, seeds, zeniths, azimuths, tel_xs, tel_zs, heights):
        total += 1
        key = key_canon(seed, energy, zen, az, h, tx, tz)
        if (not rerun_event) and key in completed:
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
echo "  RERUN_EVENT=${RERUN_EVENT}"
echo "  UPDATE_TIME_STAMP=${UPDATE_TIME_STAMP}"
if [ "${REDUCED_BASED_RUN_RADIUS_ENABLED}" = "true" ]; then
    echo "  REDUCED_BASED_RUN_RADIUS=${REDUCED_BASED_RUN_RADIUS}"
fi
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

    EXTRA_SBATCH_ARGS=()
    if [ "${REDUCED_BASED_RUN_RADIUS_ENABLED}" = "true" ]; then
        EXTRA_SBATCH_ARGS+=("--reduced_based_run_radius=${REDUCED_BASED_RUN_RADIUS}")
    fi

    SBATCH_OUTPUT=$(sbatch --job-name="corsika8_trinity_pid${PDG}_E${ENERGY}_s${SEED}" \
        "${SLURM_SCRIPT}" \
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
        "${HADRON_MODEL}" \
        "${RERUN_EVENT}" \
        --update_time_stamp="${UPDATE_TIME_STAMP}" \
        "${EXTRA_SBATCH_ARGS[@]}" 2>&1)

    if [ $? -eq 0 ]; then
        JOBID=$(echo "${SBATCH_OUTPUT}" | awk '{print $NF}')
        SUBMITTED=$((SUBMITTED + 1))

        if (( SUBMITTED % 50 == 0 )); then
            NJOBS=$(squeue -u "$USER" -h 2>/dev/null | wc -l)
            echo "  [$(date +%H:%M:%S)] Submitted: ${SUBMITTED}/${TODO_COUNT} | Queue: ${NJOBS} | rerun=${RERUN_EVENT} | E=${ENERGY} s=${SEED} zen=${ZENITH} az=${AZIMUTH} h=${INJ_HEIGHT}"
        fi
    else
        echo "    WARNING: sbatch failed: ${SBATCH_OUTPUT}"
        FAILED=$((FAILED + 1))

        if echo "${SBATCH_OUTPUT}" | grep -q "QOSMaxSubmitJobPerUserLimit"; then
            echo "    Queue limit hit. Waiting 120s and retrying..."
            sleep 120
            wait_for_queue_room

            SBATCH_OUTPUT=$(sbatch --job-name="corsika8_trinity_pid${PDG}_E${ENERGY}_s${SEED}" \
                "${SLURM_SCRIPT}" \
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
                "${HADRON_MODEL}" \
                "${RERUN_EVENT}" \
                --update_time_stamp="${UPDATE_TIME_STAMP}" \
                "${EXTRA_SBATCH_ARGS[@]}" 2>&1)

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
