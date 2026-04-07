#!/bin/bash
# =============================================================================
# submit_save_care2csv.sh
#
# Loop through multiple seeds, throttle to stay under 500 pending/running jobs,
# and produce one CSV per seed.
#
# Usage:
#   bash submit_save_care2csv.sh
# =============================================================================

CHUNK_SIZE=100
MAX_JOBS=500
SLEEP_SEC=30

PDG=13
ENERGY_MAG=5
RADIUS=5
TEL_Y=0
SEEDS=(7) #($(seq 7 26))

BASE_PATH="/scratch/general/vast/u1520754/muon_sim_chain"
ANALYSIS_DIR="$HOME/Muon_Trinity/cluster_corsika8"
WORKER_SCRIPT="${ANALYSIS_DIR}/save_CARE2csv_chunk.py"
MERGE_SCRIPT="${ANALYSIS_DIR}/merge_csv_chunks.py"

OUT_DIR="${BASE_PATH}/Muon_pid${PDG}_E1e${ENERGY_MAG}_R${RADIUS}/csv_output"

mkdir -p "$HOME/csv_logs"

for SEED in "${SEEDS[@]}"; do

    echo "============================================================"
    echo "SEED = ${SEED}"
    echo "============================================================"

    CHUNK_DIR="${OUT_DIR}/chunks_seed${SEED}"
    FINAL_CSV="${OUT_DIR}/scan_care_pid${PDG}_E1e${ENERGY_MAG}_R${RADIUS}_y${TEL_Y}_s${SEED}.csv"

    # Clear old chunks
    rm -rf "${CHUNK_DIR}"
    mkdir -p "${CHUNK_DIR}"

    # Build combos and split into chunk files
    CHUNK_SIZE_EXP=${CHUNK_SIZE} CHUNK_DIR_EXP=${CHUNK_DIR} python3 - <<'PYEOF'
import os, json
import numpy as np
from itertools import product

#zeniths  = [89.9, 89, 88, 87]
#azimuths = [268, 269, 270, 271, 272]
#tel_xs   = [-5, -2, -1, -0.7, -0.5, -0.2, 0, 0.2, 0.5, 0.7, 1, 2, 5]
#tel_zs   = [-5, -2, -1, -0.7, -0.5, -0.2, 0, 0.2, 0.5, 0.7, 1, 2, 5]
#heights  = [5000, 10000, 15000, 20000, 30000]

zeniths  = [f"{z:.1f}" for z in np.arange(87.0, 90, 0.3)]  # 87.0, 87.3, ..., 89.9
zeniths.append("89.9")  # Ensure 89.9 is included
azimuths = [f"{a:.1f}" for a in np.arange(268.0, 272.1, 0.3)]  # 268.0, 268.3, ..., 272.0
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
        echo "ERROR: Failed to create chunks for seed ${SEED}. Skipping."
        continue
    fi

    # Submit chunk array jobs in batches to stay under MAX_JOBS
    CHUNK_JOB_IDS=()
    IDX=0
    while [ ${IDX} -lt ${N_CHUNKS} ]; do

        # Wait until we have room
        while true; do
            CURRENT=$(squeue -u "$USER" -h | wc -l)
            AVAIL=$((MAX_JOBS - CURRENT))
            if [ "${AVAIL}" -gt 0 ]; then
                break
            fi
            echo "  [seed ${SEED}] ${CURRENT} jobs queued/running, waiting ${SLEEP_SEC}s ..."
            sleep ${SLEEP_SEC}
        done

        # How many chunks we can still submit
        REMAINING=$((N_CHUNKS - IDX))
        BATCH_SIZE=$((AVAIL < REMAINING ? AVAIL : REMAINING))
        END_IDX=$((IDX + BATCH_SIZE - 1))

        echo "  [seed ${SEED}] Submitting array ${IDX}-${END_IDX}  (${BATCH_SIZE} tasks)"
        mkdir -p "$HOME/csv_logs/care2csv_s${SEED}"
        JID=$(sbatch --parsable \
            --account=owner-guest \
            --partition=kingspeak-guest \
            --time=0:30:00 \
            --mem=4G \
            --array=${IDX}-${END_IDX} \
            --job-name=care2csv_s${SEED} \
            --output="$HOME/csv_logs/care2csv_s${SEED}/care2csv_s${SEED}_array_${IDX}-${END_IDX}_batch_${BATCH_SIZE}.out" \
            --error="$HOME/csv_logs/care2csv_s${SEED}/care2csv_s${SEED}_array_${IDX}-${END_IDX}_batch_${BATCH_SIZE}.err" \
            --wrap="
source /uufs/chpc.utah.edu/common/home/u1520754/miniconda3/etc/profile.d/conda.sh
conda activate jupyter_env
TASK_ID=\$(printf '%04d' \${SLURM_ARRAY_TASK_ID})
python ${WORKER_SCRIPT} \
    --chunk-file ${CHUNK_DIR}/chunk_\${TASK_ID}.json \
    --output ${CHUNK_DIR}/result_\${SLURM_ARRAY_TASK_ID}.csv \
    --pid ${PDG} \
    --energy-mag ${ENERGY_MAG} \
    --radius ${RADIUS} \
    --tel-y ${TEL_Y} \
    --seed ${SEED} \
    --base-path ${BASE_PATH}
")

        CHUNK_JOB_IDS+=("${JID}")
        IDX=$((END_IDX + 1))
    done

    # Build dependency string: afterok on all chunk array jobs for this seed
    DEP_STR=""
    for JID in "${CHUNK_JOB_IDS[@]}"; do
        if [ -z "${DEP_STR}" ]; then
            DEP_STR="afterok:${JID}"
        else
            DEP_STR="${DEP_STR}:${JID}"
        fi
    done

    # Submit merge job
    MERGE_JOB=$(sbatch --parsable \
        --account=owner-guest \
        --partition=kingspeak-guest \
        --time=0:10:00 \
        --mem=4G \
        --dependency="${DEP_STR}" \
        --job-name=merge_s${SEED} \
        --output="$HOME/csv_logs/care2csv_merge_s${SEED}.out" \
        --error="$HOME/csv_logs/care2csv_merge_s${SEED}.err" \
        --wrap="
source /uufs/chpc.utah.edu/common/home/u1520754/miniconda3/etc/profile.d/conda.sh
conda activate jupyter_env
python3 ${MERGE_SCRIPT} \
    --chunk-dir ${CHUNK_DIR} \
    --output ${FINAL_CSV}
")

    echo "  Merge job: ${MERGE_JOB} (depends on ${DEP_STR})"
    echo "  Final CSV: ${FINAL_CSV}"
    echo ""

done

echo "All seeds submitted."
