#!/bin/bash
CHUNK_SIZE=100
MAX_JOBS=500
SLEEP_SEC=30
MAX_SUBMIT_RETRIES=20

PDG=""
RADIUS=""
TEL_Y=""
SEEDS=()
ONLY_ENERGY=""
ONLY_SEED=""
WAIT_FOR_MERGE="false"
MAX_PE_CUT=""
DRY_RUN_DELETE="false"

BASE_PATH="/scratch/general/vast/u1520754/muon_sim_chain_tree"
ANALYSIS_DIR="$HOME/Muon_Trinity/cluster_corsika8/save_data2csv"
WORKER_SCRIPT="${ANALYSIS_DIR}/save_CARE2csv_chunk_tree.py"
MERGE_SCRIPT="${ANALYSIS_DIR}/merge_csv_chunks.py"
CORRECTION_REPORT_NAME="metadata.yaml"
SIM_INPUT_SCRIPT="$HOME/Muon_Trinity/cluster_corsika8/Trinity_sim_imput.py"

while [[ $# -gt 0 ]]; do
    case "$1" in
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
        --wait|--wait-merge)
            WAIT_FOR_MERGE="true"
            shift
            ;;
        --Max_PE_cut)
            if [[ $# -gt 1 && ! "$2" =~ ^-- ]]; then
                MAX_PE_CUT="$2"
                shift 2
            else
                echo "ERROR: --Max_PE_cut requires a value"
                exit 1
            fi
            ;;
        --Max_PE_cut=*)
            MAX_PE_CUT="${1#*=}"
            shift
            ;;
        --dry-run-delete)
            DRY_RUN_DELETE="true"
            shift
            ;;
        --help|-h)
            echo "Usage: bash submit_save_care2csv_tree.sh [--only-energy E] [--only-seed S] [--wait] [--Max_PE_cut V] [--dry-run-delete]"
            echo "  --only-energy E   Process only a single energy string"
            echo "  --only-seed S     Process only a single seed"
            echo "  --wait            Block until merge job finishes"
            echo "  --Max_PE_cut V    Delete run dir when max_pe < V (PE)"
            echo "  --dry-run-delete  Log low-PE deletions without removing files"
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option '$1'"
            echo "Usage: bash submit_save_care2csv_tree.sh [--only-energy E] [--only-seed S] [--wait]"
            exit 1
            ;;
    esac
done

mkdir -p "$HOME/csv_logs"

source /uufs/chpc.utah.edu/common/home/u1520754/miniconda3/etc/profile.d/conda.sh
conda activate jupyter_env
SIM_INPUT_JSON=$(mktemp /tmp/trinity_sim_input_XXXXXX.json)
python3 "${SIM_INPUT_SCRIPT}" --output "${SIM_INPUT_JSON}"

read -r PDG TEL_Y RADIUS SEEDS_CSV <<< "$(python3 - "${SIM_INPUT_JSON}" <<'PYEOF'
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
python3 - "${SIM_INPUT_JSON}" <<'PYEOF' > "$ENERGY_FILE"
import json
import sys

with open(sys.argv[1], "r", encoding="utf-8") as f:
    data = json.load(f)

for energy in data.get("energy_strings", []):
    print(energy)
PYEOF

mapfile -t ENERGY_STRS < "$ENERGY_FILE"
rm -f "$ENERGY_FILE"

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
    ENERGY_STRS=("${ONLY_ENERGY}")
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
    SEEDS=("${ONLY_SEED}")
fi

echo "Found ${#ENERGY_STRS[@]} energies from MCEq grid"
echo "Energies: ${ENERGY_STRS[*]}"
echo "CSV completion rule: requires BOTH CARE/cherenkov_hits.root and ${CORRECTION_REPORT_NAME}"
if [ -n "${MAX_PE_CUT}" ]; then
    echo "  Low-PE deletion: max_pe < ${MAX_PE_CUT} marks file_found=1 and deletes run dir"
    if [ "${DRY_RUN_DELETE}" = "true" ]; then
        echo "  DRY-RUN enabled: no deletions will occur"
    fi
fi

wait_for_job_completion() {
    local job_id="$1"
    if [ -z "${job_id}" ]; then
        return 1
    fi
    while squeue -j "${job_id}" -h | grep -q .; do
        echo "  [merge ${job_id}] still running, waiting ${SLEEP_SEC}s ..."
        sleep ${SLEEP_SEC}
    done
    return 0
}

# --- Helper: submit with retry on QOS/transient errors ---
submit_with_retry() {
    local retries=0
    local JID=""
    while [ $retries -lt $MAX_SUBMIT_RETRIES ]; do
        JID=$(sbatch --parsable "$@" 2>/tmp/sbatch_err_$$)
        local rc=$?
        local errmsg=$(cat /tmp/sbatch_err_$$ 2>/dev/null)
        rm -f /tmp/sbatch_err_$$

        if [ $rc -eq 0 ] && [ -n "$JID" ]; then
            echo "$JID"
            return 0
        fi

        # Check if it's a QOS limit or transient error worth retrying
        if echo "$errmsg" | grep -qi "QOSMaxSubmitJobPerUserLimit\|Resource temporarily unavailable\|Socket timed out"; then
            retries=$((retries + 1))
            echo "  [RETRY $retries/$MAX_SUBMIT_RETRIES] QOS/transient limit hit, waiting ${SLEEP_SEC}s ..." >&2
            sleep ${SLEEP_SEC}
        else
            # Non-retryable error
            echo "  [ERROR] sbatch failed: $errmsg" >&2
            echo ""
            return 1
        fi
    done
    echo "  [ERROR] Exhausted $MAX_SUBMIT_RETRIES retries" >&2
    echo ""
    return 1
}

# --- Count total submitted array elements (not just job lines) ---
count_submitted_elements() {
    # Each array job line from squeue may represent many tasks
    # Use squeue to count individual array tasks
    squeue -u "$USER" -h -r | wc -l
}

for ENERGY_STR in "${ENERGY_STRS[@]}"; do
for SEED in "${SEEDS[@]}"; do

    echo "============================================================"
    echo "ENERGY = ${ENERGY_STR}, SEED = ${SEED}"
    echo "============================================================"

    OUT_DIR="${BASE_PATH}/Muon_pid${PDG}_E${ENERGY_STR}_R${RADIUS}/csv_output"
    CHUNK_DIR="${OUT_DIR}/chunks_seed${SEED}"
    FINAL_CSV="${OUT_DIR}/scan_care_pid${PDG}_E${ENERGY_STR}_R${RADIUS}_y${TEL_Y}_s${SEED}.csv"

    rm -rf "${CHUNK_DIR}"
    mkdir -p "${CHUNK_DIR}"

SIM_INPUT_JSON="${SIM_INPUT_JSON}" CHUNK_SIZE_EXP=${CHUNK_SIZE} CHUNK_DIR_EXP=${CHUNK_DIR} python3 - <<'PYEOF'
import os, json
from itertools import product

sim_input_path = os.environ.get("SIM_INPUT_JSON")
if not sim_input_path:
    raise SystemExit("SIM_INPUT_JSON is not set")

with open(sim_input_path, "r", encoding="utf-8") as f:
    sim = json.load(f)

geom = sim.get("geometry", {})
zeniths = geom.get("zeniths_deg", [])
azimuths = geom.get("azimuths_deg", [])
tel_xs = geom.get("tel_xs_m", [])
tel_zs = geom.get("tel_zs_m", [])
heights = geom.get("heights_m", [])

combos = list(product(zeniths, azimuths, heights, tel_xs, tel_zs))
chunk_size = int(os.environ["CHUNK_SIZE_EXP"])
chunk_dir  = os.environ["CHUNK_DIR_EXP"]

n_chunks = 0
for i in range(0, len(combos), chunk_size):
    chunk = combos[i:i+chunk_size]
    chunk_file = os.path.join(chunk_dir, f"chunk_{n_chunks:04d}.json")
    with open(chunk_file, "w") as f:
        json.dump(chunk, f)
    n_chunks += 1

print(f"  Total combinations: {len(combos)}")
print(f"  Chunks created: {n_chunks}")

with open(os.path.join(chunk_dir, "n_chunks.txt"), "w") as f:
    f.write(str(n_chunks))
PYEOF

    N_CHUNKS=$(cat "${CHUNK_DIR}/n_chunks.txt")
    if [ -z "${N_CHUNKS}" ] || [ "${N_CHUNKS}" -lt 1 ]; then
        echo "ERROR: Failed to create chunks. Skipping."
        continue
    fi

    CHUNK_JOB_IDS=()
    ALL_SUBMITTED=true
    IDX=0
    while [ ${IDX} -lt ${N_CHUNKS} ]; do
        # Wait until we have room
        while true; do
            CURRENT=$(count_submitted_elements)
            AVAIL=$((MAX_JOBS - CURRENT))
            [ "${AVAIL}" -gt 10 ] && break  # need at least some headroom
            echo "  [E=${ENERGY_STR} s${SEED}] ${CURRENT} elements queued, waiting ${SLEEP_SEC}s ..."
            sleep ${SLEEP_SEC}
        done

        REMAINING=$((N_CHUNKS - IDX))
        # Be conservative: leave buffer of 10 for merge jobs etc.
        SAFE_AVAIL=$((AVAIL - 10))
        [ "${SAFE_AVAIL}" -lt 1 ] && SAFE_AVAIL=1
        BATCH_SIZE=$((SAFE_AVAIL < REMAINING ? SAFE_AVAIL : REMAINING))
        END_IDX=$((IDX + BATCH_SIZE - 1))

        echo "  [E=${ENERGY_STR} s${SEED}] Submitting array ${IDX}-${END_IDX} (${BATCH_SIZE} tasks)"
        LOG_DIR="$HOME/csv_logs/care2csv_E${ENERGY_STR}_s${SEED}"
        mkdir -p "$LOG_DIR"

        EXTRA_WORKER_ARGS=""
        if [ -n "${MAX_PE_CUT}" ]; then
            EXTRA_WORKER_ARGS+=" --Max_PE_cut ${MAX_PE_CUT}"
        fi
        if [ "${DRY_RUN_DELETE}" = "true" ]; then
            EXTRA_WORKER_ARGS+=" --dry-run-delete"
        fi

        JID=$(submit_with_retry \
            --account=owner-guest \
            --partition=kingspeak-guest \
            --time=0:30:00 \
            --mem=4G \
            --array=${IDX}-${END_IDX} \
            --job-name=care2csv_E${ENERGY_STR}_s${SEED} \
            --output="${LOG_DIR}/array_%a.out" \
            --error="${LOG_DIR}/array_%a.err" \
            --wrap="
source /uufs/chpc.utah.edu/common/home/u1520754/miniconda3/etc/profile.d/conda.sh
conda activate jupyter_env
TASK_ID=\$(printf '%04d' \${SLURM_ARRAY_TASK_ID})
python ${WORKER_SCRIPT} \
    --chunk-file ${CHUNK_DIR}/chunk_\${TASK_ID}.json \
    --output ${CHUNK_DIR}/result_\${SLURM_ARRAY_TASK_ID}.csv \
    --pid ${PDG} \
    --energy-str ${ENERGY_STR} \
    --radius ${RADIUS} \
    --tel-y ${TEL_Y} \
    --seed ${SEED} \
    --correction-report-name ${CORRECTION_REPORT_NAME} \
    --base-path ${BASE_PATH}${EXTRA_WORKER_ARGS}
")

        if [ -n "$JID" ]; then
            CHUNK_JOB_IDS+=("${JID}")
            IDX=$((END_IDX + 1))
        else
            echo "  [E=${ENERGY_STR} s${SEED}] ERROR: Failed to submit array ${IDX}-${END_IDX} after retries"
            ALL_SUBMITTED=false
            break
        fi

        # Small delay to let scheduler update
        sleep 2
    done

    if [ "$ALL_SUBMITTED" = false ] || [ ${#CHUNK_JOB_IDS[@]} -eq 0 ]; then
        echo "  [E=${ENERGY_STR} s${SEED}] SKIPPING merge — not all chunks submitted."
        echo ""
        continue
    fi

    DEP_STR="afterany"
    for JID in "${CHUNK_JOB_IDS[@]}"; do
        DEP_STR="${DEP_STR}:${JID}"
    done

    MERGE_JOB=$(submit_with_retry \
        --account=owner-guest \
        --partition=kingspeak-guest \
        --time=0:10:00 \
        --mem=4G \
        --dependency="${DEP_STR}" \
        --job-name=merge_E${ENERGY_STR}_s${SEED} \
        --output="$HOME/csv_logs/merge_E${ENERGY_STR}_s${SEED}.out" \
        --error="$HOME/csv_logs/merge_E${ENERGY_STR}_s${SEED}.err" \
        --wrap="
source /uufs/chpc.utah.edu/common/home/u1520754/miniconda3/etc/profile.d/conda.sh
conda activate jupyter_env
python3 ${MERGE_SCRIPT} \
    --chunk-dir ${CHUNK_DIR} \
    --output ${FINAL_CSV}
")

    if [ -n "$MERGE_JOB" ]; then
        echo "  Merge job: ${MERGE_JOB} -> ${FINAL_CSV}"
        if [ "${WAIT_FOR_MERGE}" = "true" ]; then
            wait_for_job_completion "${MERGE_JOB}"
            echo "  Merge completed: ${FINAL_CSV}"
        fi
    else
        echo "  [E=${ENERGY_STR} s${SEED}] ERROR: Failed to submit merge job"
    fi
    echo ""

done
done

echo "All seeds/energies submitted."
