#!/bin/bash
# submit_trigger_resample.sh — chunk triggered events and submit SLURM array jobs

CHUNK_SIZE=100
MAX_JOBS=500
SLEEP_SEC=30
MAX_SUBMIT_RETRIES=20

PDG=13
RADIUS=5
TEL_Y=0
PE_THRESHOLD=20
SEEDS=(1 2 3)

TEL_X_VALUES="-5 -2 -1 -0.3 0 0.3 1 2 5"
TEL_Z_VALUES="-5 -2 -1 -0.3 0 0.3 1 2 5"

BASE_PATH="/scratch/general/vast/u1520754/muon_sim_chain_tree"
ANALYSIS_DIR="$HOME/Muon_Trinity/cluster_corsika8/save_data2csv"
WORKER_SCRIPT="${ANALYSIS_DIR}/build_trigger_resample_chunk.py"
MERGE_SCRIPT="${ANALYSIS_DIR}/merge_trigger_resample_csv.py"

mkdir -p "$HOME/csv_logs"

# --- Generate energy strings from MCEq ---
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
E = E[(E >= 1e3) & (E <= 4e3)]

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

# --- Helper: submit with retry ---
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

        if echo "$errmsg" | grep -qi "QOSMaxSubmitJobPerUserLimit\|Resource temporarily unavailable\|Socket timed out"; then
            retries=$((retries + 1))
            echo "  [RETRY $retries/$MAX_SUBMIT_RETRIES] waiting ${SLEEP_SEC}s ..." >&2
            sleep ${SLEEP_SEC}
        else
            echo "  [ERROR] sbatch failed: $errmsg" >&2
            echo ""
            return 1
        fi
    done
    echo "  [ERROR] Exhausted $MAX_SUBMIT_RETRIES retries" >&2
    echo ""
    return 1
}

# --- Main loop ---
for SEED in "${SEEDS[@]}"; do
for ENERGY_STR in "${ENERGY_STRS[@]}"; do

    ENERGY_GROUP="Muon_pid${PDG}_E${ENERGY_STR}_R${RADIUS}"
    CSV_DIR="${BASE_PATH}/${ENERGY_GROUP}/csv_output"
    ORIG_CSV="${CSV_DIR}/scan_care_pid${PDG}_E${ENERGY_STR}_R${RADIUS}_y${TEL_Y}_s${SEED}.csv"
    FINAL_CSV="${CSV_DIR}/scan_care_pid${PDG}_E${ENERGY_STR}_R${RADIUS}_y${TEL_Y}_s${SEED}_trigger.csv"

    if [ ! -f "$ORIG_CSV" ]; then
        echo "[E=${ENERGY_STR} s${SEED}] Original CSV not found, skipping: ${ORIG_CSV}"
        continue
    fi

    if [ -f "$FINAL_CSV" ]; then
        echo "[E=${ENERGY_STR} s${SEED}] Output exists, will overwrite: ${FINAL_CSV}"
        rm -f "$FINAL_CSV"
    fi


    echo "[E=${ENERGY_STR} s${SEED}] Filtering triggered events (PE >= ${PE_THRESHOLD})..."

    # --- Step 1: Filter triggered events and write chunk JSON files ---
    CHUNK_DIR="${CSV_DIR}/trigger_chunks_s${SEED}"
    rm -rf "$CHUNK_DIR"
    mkdir -p "$CHUNK_DIR"

    # Python: read CSV, filter, write chunks
    N_CHUNKS=$(python3 - "$ORIG_CSV" "$CHUNK_DIR" "$PE_THRESHOLD" "$CHUNK_SIZE" <<'FILTER_EOF'
import sys, csv, json, os, math

orig_csv = sys.argv[1]
chunk_dir = sys.argv[2]
pe_threshold = float(sys.argv[3])
chunk_size = int(sys.argv[4])

triggered = []
with open(orig_csv, "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            if float(row["max_pe"]) >= pe_threshold:
                triggered.append(row)
        except (ValueError, KeyError):
            pass

n_chunks = max(1, math.ceil(len(triggered) / chunk_size))
for i in range(n_chunks):
    chunk = triggered[i*chunk_size : (i+1)*chunk_size]
    with open(os.path.join(chunk_dir, f"chunk_{i:04d}.json"), "w") as jf:
        json.dump(chunk, jf)

print(n_chunks)
FILTER_EOF
    )

    if [ -z "$N_CHUNKS" ] || [ "$N_CHUNKS" -eq 0 ]; then
        echo "  [E=${ENERGY_STR} s${SEED}] No triggered events found, skipping."
        continue
    fi

    echo "  [E=${ENERGY_STR} s${SEED}] ${N_CHUNKS} chunks created"

    # --- Step 2: Submit SLURM array jobs in batches ---
    CHUNK_JOB_IDS=()
    IDX=0
    MAX_IDX=$((N_CHUNKS - 1))
    ALL_SUBMITTED=true

    while [ $IDX -le $MAX_IDX ]; do
        END_IDX=$((IDX + MAX_JOBS - 1))
        if [ $END_IDX -gt $MAX_IDX ]; then
            END_IDX=$MAX_IDX
        fi

        BATCH_SIZE=$((END_IDX - IDX + 1))
        echo "  [E=${ENERGY_STR} s${SEED}] Submitting array ${IDX}-${END_IDX} (${BATCH_SIZE} tasks)"
        LOG_DIR="$HOME/csv_logs/trigger_E${ENERGY_STR}_s${SEED}"
        mkdir -p "$LOG_DIR"

        JID=$(submit_with_retry \
            --account=owner-guest \
            --partition=kingspeak-guest \
            --time=1:00:00 \
            --mem=4G \
            --array=${IDX}-${END_IDX} \
            --job-name=trig_E${ENERGY_STR}_s${SEED} \
            --output="${LOG_DIR}/array_%a.out" \
            --error="${LOG_DIR}/array_%a.err" \
            --wrap="
source /uufs/chpc.utah.edu/common/home/u1520754/miniconda3/etc/profile.d/conda.sh
conda activate jupyter_env
TASK_ID=\$(printf '%04d' \${SLURM_ARRAY_TASK_ID})
python3 ${WORKER_SCRIPT} \
    --chunk-file ${CHUNK_DIR}/chunk_\${TASK_ID}.json \
    --output ${CHUNK_DIR}/result_\${SLURM_ARRAY_TASK_ID}.csv \
    --pid ${PDG} \
    --energy-str ${ENERGY_STR} \
    --radius ${RADIUS} \
    --tel-y ${TEL_Y} \
    --seed ${SEED} \
    --base-path ${BASE_PATH} \
    --tel-x-values ${TEL_X_VALUES} \
    --tel-z-values ${TEL_Z_VALUES}
")

        if [ -n "$JID" ]; then
            CHUNK_JOB_IDS+=("${JID}")
            IDX=$((END_IDX + 1))
        else
            echo "  [E=${ENERGY_STR} s${SEED}] ERROR: Failed to submit array ${IDX}-${END_IDX}"
            ALL_SUBMITTED=false
            break
        fi

        sleep 2
    done

    if [ "$ALL_SUBMITTED" = false ] || [ ${#CHUNK_JOB_IDS[@]} -eq 0 ]; then
        echo "  [E=${ENERGY_STR} s${SEED}] SKIPPING merge — not all chunks submitted."
        continue
    fi

    # --- Step 3: Submit merge job ---
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
        --job-name=merge_trig_E${ENERGY_STR}_s${SEED} \
        --output="$HOME/csv_logs/merge_trig_E${ENERGY_STR}_s${SEED}.out" \
        --error="$HOME/csv_logs/merge_trig_E${ENERGY_STR}_s${SEED}.err" \
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
