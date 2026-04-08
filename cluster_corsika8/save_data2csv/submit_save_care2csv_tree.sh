#!/bin/bash
CHUNK_SIZE=100
MAX_JOBS=500
SLEEP_SEC=30
MAX_SUBMIT_RETRIES=20

PDG=13
RADIUS=5
TEL_Y=0
SEEDS=(1)

BASE_PATH="/scratch/general/vast/u1520754/muon_sim_chain_tree"
ANALYSIS_DIR="$HOME/Muon_Trinity/cluster_corsika8/save_data2csv"
WORKER_SCRIPT="${ANALYSIS_DIR}/save_CARE2csv_chunk_tree.py"
MERGE_SCRIPT="${ANALYSIS_DIR}/merge_csv_chunks.py"

mkdir -p "$HOME/csv_logs"

# Generate energy strings from MCEq
source /uufs/chpc.utah.edu/common/home/u1520754/miniconda3/etc/profile.d/conda.sh
conda activate jupyter_env

ENERGY_FILE=$(mktemp)
python3 <<'PYEOF' > "$ENERGY_FILE"
import sys, os, math
os.environ["MCEQ_LOG_LEVEL"] = "50"
import logging
logging.disable(logging.CRITICAL)
import MCEq.config as config
config.kernel_config = 'MKL'
config.e_min = 1000
config.integrator = 'euler'

old_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')
from MCEq.core import MCEqRun
import crflux.models as pm
mceq = MCEqRun(
    interaction_model='SIBYLL23C',
    primary_model=(pm.GlobalSplineFitBeta, None),
    theta_deg=0.,
    density_model=("CORSIKA", ('USStd', None)),
)
sys.stdout = old_stdout

E = mceq.e_grid
E = E[(E > 5e3) & (E <= 1e4)]

def to_1e(val):
    exp = int(math.log10(val))
    coeff = val / 10**exp
    return f"{coeff:g}e{exp}"

for e in E:
    print(to_1e(e))
PYEOF

mapfile -t ENERGY_STRS < "$ENERGY_FILE"
rm -f "$ENERGY_FILE"

echo "Found ${#ENERGY_STRS[@]} energies from MCEq grid"
echo "Energies: ${ENERGY_STRS[*]}"

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

    CHUNK_SIZE_EXP=${CHUNK_SIZE} CHUNK_DIR_EXP=${CHUNK_DIR} python3 - <<'PYEOF'
import os, json
import numpy as np
from itertools import product

zeniths  = [f"{z:.1f}" for z in np.arange(87.0, 90, 0.3)]
zeniths.append("89.9")
azimuths = [f"{a:.1f}" for a in np.arange(267.0, 273.3, 0.3)]
tel_xs   = ["0"]
tel_zs   = ["0"]
heights  = [str(int(h)) for h in np.linspace(5000, 50000, 100)]

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
    --base-path ${BASE_PATH}
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
    else
        echo "  [E=${ENERGY_STR} s${SEED}] ERROR: Failed to submit merge job"
    fi
    echo ""

done
done

echo "All seeds/energies submitted."
