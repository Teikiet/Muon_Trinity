#!/bin/bash
# =============================================================================
# submit_save_care2csv.sh
# =============================================================================

CHUNK_SIZE=500

PDG=13
ENERGY_MAG=5
RADIUS=5
TEL_Y=0
SEED=4

BASE_PATH="/scratch/general/vast/u1520754/muon_sim_chain"
ANALYSIS_DIR="$HOME/Muon_Trinity/cluster_corsika8"
WORKER_SCRIPT="${ANALYSIS_DIR}/save_CARE2csv_chunk.py"
MERGE_SCRIPT="${ANALYSIS_DIR}/merge_csv_chunks.py"

OUT_DIR="${BASE_PATH}/Muon_pid${PDG}_E1e${ENERGY_MAG}_R${RADIUS}/csv_output"
CHUNK_DIR="${OUT_DIR}/chunks_seed${SEED}"
# Clear old chunks before starting
rm -rf "${CHUNK_DIR}"
mkdir -p "${CHUNK_DIR}"

FINAL_CSV="${OUT_DIR}/scan_care_pid${PDG}_E1e${ENERGY_MAG}_R${RADIUS}_y${TEL_Y}_s${SEED}.csv"

# Export so the heredoc Python can see them
export CHUNK_SIZE CHUNK_DIR

# Build all combos and split into chunk files
python - <<'PYEOF'
import os, json
from itertools import product

zeniths  = [89.9, 89, 88, 87]
azimuths = [268, 269, 270, 271, 272]
tel_xs   = [-5, -2, -1, -0.7, -0.5, -0.2, 0, 0.2, 0.5, 0.7, 1, 2, 5]
tel_zs   = [-5, -2, -1, -0.7, -0.5, -0.2, 0, 0.2, 0.5, 0.7, 1, 2, 5]
heights  = [5000, 10000, 15000, 20000, 30000]

combos = list(product(zeniths, azimuths, heights, tel_xs, tel_zs))
chunk_size = int(os.environ["CHUNK_SIZE"])
chunk_dir  = os.environ["CHUNK_DIR"]

n_chunks = 0
for i in range(0, len(combos), chunk_size):
    chunk = combos[i:i+chunk_size]
    chunk_file = os.path.join(chunk_dir, f"chunk_{n_chunks:04d}.json")
    with open(chunk_file, "w") as f:
        json.dump(chunk, f)
    n_chunks += 1

print(f"Total combinations: {len(combos)}")
print(f"Chunks created: {n_chunks}")

with open(os.path.join(chunk_dir, "n_chunks.txt"), "w") as f:
    f.write(str(n_chunks))
PYEOF

N_CHUNKS=$(cat "${CHUNK_DIR}/n_chunks.txt")

if [ -z "${N_CHUNKS}" ] || [ "${N_CHUNKS}" -lt 1 ]; then
    echo "ERROR: Failed to create chunks. Exiting."
    exit 1
fi

echo "============================================================"
echo "Submitting ${N_CHUNKS} chunk jobs + 1 merge job"
echo "Output: ${FINAL_CSV}"
echo "============================================================"

# Pad array task ID to 4 digits in the --wrap command
ARRAY_JOB=$(sbatch --parsable \
    --account=owner-guest \
    --partition=kingspeak-guest \
    --time=0:30:00 \
    --mem=4G \
    --array=0-$((N_CHUNKS - 1)) \
    --job-name=care2csv_chunk \
    --output="$HOME/log/care2csv_%A_%a.out" \
    --error="$HOME/log/care2csv_%A_%a.err" \
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

echo "Chunk array job: ${ARRAY_JOB}"

# Submit merge job with dependency
MERGE_JOB=$(sbatch --parsable \
    --account=owner-guest \
    --partition=kingspeak-guest \
    --time=0:10:00 \
    --mem=4G \
    --dependency=afterok:${ARRAY_JOB} \
    --job-name=care2csv_merge \
    --output="$HOME/log/care2csv_merge_%j.out" \
    --error="$HOME/log/care2csv_merge_%j.err" \
    --wrap="
source /uufs/chpc.utah.edu/common/home/u1520754/miniconda3/etc/profile.d/conda.sh
conda activate jupyter_env
python3 ${MERGE_SCRIPT} \
    --chunk-dir ${CHUNK_DIR} \
    --output ${FINAL_CSV}
")

echo "Merge job: ${MERGE_JOB} (depends on ${ARRAY_JOB})"
echo "Final CSV: ${FINAL_CSV}"
