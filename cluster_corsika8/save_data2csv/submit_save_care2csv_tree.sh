#!/bin/bash
CHUNK_SIZE=100
MAX_JOBS=500
SLEEP_SEC=30
MAX_SUBMIT_RETRIES=20
UPLOAD_TRIGGERED_TO_DRIVE=""
UPLOAD_TRIGGERED_THRESHOLD=20

PDG=""
RADIUS=""
TEL_Y=""
SEEDS=()
ONLY_ENERGY=""
ONLY_SEED=""
WAIT_FOR_MERGE="false"
MAX_PE_CUT=""
DRY_RUN_DELETE="false"
TRIGGERED_BASE_ONLY="false"
TRIGGERED_BASE_MAX_PE="20"

BASE_PATH="/scratch/general/vast/u1520754/muon_sim_chain_tree"
ANALYSIS_DIR="$HOME/Muon_Trinity/cluster_corsika8/save_data2csv"
WORKER_SCRIPT="${ANALYSIS_DIR}/save_CARE2csv_chunk_tree.py"
MERGE_SCRIPT="${ANALYSIS_DIR}/merge_csv_chunks.py"
UPLOAD_SCRIPT="${ANALYSIS_DIR}/upload_triggered_care_to_drive.py"
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
        --triggered-base)
            TRIGGERED_BASE_ONLY="true"
            shift
            ;;
        --triggered-base-only)
            TRIGGERED_BASE_ONLY="true"
            shift
            ;;
        --triggered-base-max-pe)
            if [[ $# -gt 1 && ! "$2" =~ ^-- ]]; then
                TRIGGERED_BASE_MAX_PE="$2"
                shift 2
            else
                echo "ERROR: --triggered-base-max-pe requires a value"
                exit 1
            fi
            ;;
        --triggered-base-max-pe=*)
            TRIGGERED_BASE_MAX_PE="${1#*=}"
            shift
            ;;
        --upload-triggered-drive)
            if [[ $# -gt 1 && ! "$2" =~ ^-- ]]; then
                UPLOAD_TRIGGERED_TO_DRIVE="$2"
                shift 2
            else
                echo "ERROR: --upload-triggered-drive requires a value"
                exit 1
            fi
            ;;
        --upload-triggered-drive=*)
            UPLOAD_TRIGGERED_TO_DRIVE="${1#*=}"
            shift
            ;;
        --upload-triggered-threshold)
            if [[ $# -gt 1 && ! "$2" =~ ^-- ]]; then
                UPLOAD_TRIGGERED_THRESHOLD="$2"
                shift 2
            else
                echo "ERROR: --upload-triggered-threshold requires a value"
                exit 1
            fi
            ;;
        --upload-triggered-threshold=*)
            UPLOAD_TRIGGERED_THRESHOLD="${1#*=}"
            shift
            ;;
        --help|-h)
            echo "Usage: bash submit_save_care2csv_tree.sh [--only-energy E] [--only-seed S] [--wait] [--Max_PE_cut V] [--dry-run-delete] [--triggered-base] [--triggered-base-max-pe PE] [--upload-triggered-drive DEST] [--upload-triggered-threshold PE]"
            echo "  --only-energy E   Process only a single energy string"
            echo "  --only-seed S     Process only a single seed"
            echo "  --wait            Block until merge job finishes"
            echo "  --Max_PE_cut V    Delete run dir when max_pe < V (PE)"
            echo "  --dry-run-delete  Log low-PE deletions without removing files"
            echo "  --triggered-base   Only process offset geometries whose base row has max_pe >= threshold"
            echo "  --triggered-base-max-pe PE   Trigger threshold for base rows (default: 20)"
            echo "  --upload-triggered-drive DEST   Upload triggered CARE roots to an rclone destination"
            echo "  --upload-triggered-threshold PE  Trigger threshold for uploads (default: 20)"
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option '$1'"
            echo "Usage: bash submit_save_care2csv_tree.sh [--only-energy E] [--only-seed S] [--wait] [--triggered-base]"
            exit 1
            ;;
    esac
done

if ! [[ "${TRIGGERED_BASE_MAX_PE}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "ERROR: --triggered-base-max-pe must be a positive number. Got '${TRIGGERED_BASE_MAX_PE}'"
    exit 1
fi

if ! awk "BEGIN {exit !(${TRIGGERED_BASE_MAX_PE} > 0)}"; then
    echo "ERROR: --triggered-base-max-pe must be > 0. Got '${TRIGGERED_BASE_MAX_PE}'"
    exit 1
fi

if [ -n "${UPLOAD_TRIGGERED_TO_DRIVE}" ]; then
    if ! command -v rclone >/dev/null 2>&1; then
        echo "ERROR: --upload-triggered-drive requires rclone to be available"
        exit 1
    fi
    if ! awk "BEGIN {exit !(${UPLOAD_TRIGGERED_THRESHOLD} > 0)}"; then
        echo "ERROR: --upload-triggered-threshold must be > 0. Got '${UPLOAD_TRIGGERED_THRESHOLD}'"
        exit 1
    fi
fi

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
if [ "${TRIGGERED_BASE_ONLY}" = "true" ]; then
    echo "Triggered-base mode: only offset geometries with base max_pe >= ${TRIGGERED_BASE_MAX_PE} will be saved"
fi
if [ -n "${MAX_PE_CUT}" ]; then
    echo "  Low-PE deletion: max_pe < ${MAX_PE_CUT} marks file_found=1 and deletes run dir"
    if [ "${DRY_RUN_DELETE}" = "true" ]; then
        echo "  DRY-RUN enabled: no deletions will occur"
    fi
fi
if [ -n "${UPLOAD_TRIGGERED_TO_DRIVE}" ]; then
    echo "Triggered upload: max_pe >= ${UPLOAD_TRIGGERED_THRESHOLD} will be copied to ${UPLOAD_TRIGGERED_TO_DRIVE}"
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

build_triggered_base_chunk_manifest() {
    local energy_str="$1"
    local seed="$2"
    local manifest_file="$3"

    python3 - "$SIM_INPUT_JSON" "$BASE_PATH" "$PDG" "$RADIUS" "$TEL_Y" "$energy_str" "$seed" "$TRIGGERED_BASE_MAX_PE" "$manifest_file" <<'PYEOF'
import csv
import json
import os
import sys

sim_input_path, base_path, pdg, radius, tel_y, energy_str, seed, max_pe_cut, manifest_file = sys.argv[1:]
cut = float(max_pe_cut)

with open(sim_input_path, "r", encoding="utf-8") as f:
    sim = json.load(f)

geom = sim.get("geometry", {})
zeniths = geom.get("zeniths_deg", [])
azimuths = geom.get("azimuths_deg", [])
tel_xs = geom.get("tel_xs_m", [])
tel_zs = geom.get("tel_zs_m", [])
heights = geom.get("heights_m", [])

csv_path = os.path.join(
    base_path,
    f"Muon_pid{pdg}_E{energy_str}_R{radius}",
    "csv_output",
    f"scan_care_pid{pdg}_E{energy_str}_R{radius}_y{tel_y}_s{seed}.csv",
)

triggered = set()
if os.path.exists(csv_path):
    with open(csv_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                if float(row.get("tel_x", "nan")) != 0.0:
                    continue
                if float(row.get("tel_z", "nan")) != 0.0:
                    continue
                if str(row.get("file_found", "0")).strip() != "1":
                    continue
                if float(row.get("max_pe", "nan")) < cut:
                    continue
                triggered.add((
                    str(int(float(row["seed"]))),
                    f"{float(row['zen']):.1f}",
                    f"{float(row['az']):.1f}",
                    str(int(round(float(row['height'])))),
                ))
            except Exception:
                continue

total = 0
selected = 0
with open(manifest_file, "w", encoding="utf-8") as out:
    for zen, az, height, tel_x, tel_z in __import__("itertools").product(zeniths, azimuths, heights, tel_xs, tel_zs):
        total += 1
        if float(tel_x) == 0.0 and float(tel_z) == 0.0:
            continue
        base_key = (
            str(int(float(seed))),
            f"{float(zen):.1f}",
            f"{float(az):.1f}",
            str(int(round(float(height)))),
        )
        if base_key not in triggered:
            continue
        out.write(f"{zen}\t{az}\t{height}\t{tel_x}\t{tel_z}\n")
        selected += 1

print(f"triggered={len(triggered)} selected={selected} total={total}")
PYEOF
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

    if [ "${TRIGGERED_BASE_ONLY}" = "true" ]; then
        MANIFEST_FILE="${CHUNK_DIR}/triggered_manifest.tsv"
        build_triggered_base_chunk_manifest "${ENERGY_STR}" "${SEED}" "${MANIFEST_FILE}"
        if [ ! -s "${MANIFEST_FILE}" ]; then
            echo "  No triggered base rows found for E=${ENERGY_STR} s=${SEED}; skipping CSV jobs."
            rm -rf "${CHUNK_DIR}"
            continue
        fi
    fi

SIM_INPUT_JSON="${SIM_INPUT_JSON}" CHUNK_SIZE_EXP=${CHUNK_SIZE} CHUNK_DIR_EXP=${CHUNK_DIR} TRIGGERED_BASE_ONLY="${TRIGGERED_BASE_ONLY}" MANIFEST_FILE_EXP="${CHUNK_DIR}/triggered_manifest.tsv" python3 - <<'PYEOF'
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
triggered_base_only = os.environ.get("TRIGGERED_BASE_ONLY", "false").strip().lower() in ("true", "1", "yes", "y")
manifest_file = os.environ.get("MANIFEST_FILE_EXP", "")
allowed = None
if triggered_base_only and manifest_file and os.path.exists(manifest_file):
    allowed = set()
    with open(manifest_file, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != 5:
                continue
            allowed.add(tuple(parts))

if allowed is not None:
    filtered = []
    for zen, az, height, tel_x, tel_z in combos:
        if float(tel_x) == 0.0 and float(tel_z) == 0.0:
            continue
        key = (
            f"{float(zen):.1f}",
            f"{float(az):.1f}",
            str(int(round(float(height)))),
            str(tel_x),
            str(tel_z),
        )
        if key in allowed:
            filtered.append((zen, az, height, tel_x, tel_z))
    combos = filtered

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
        if [ -n "${UPLOAD_TRIGGERED_TO_DRIVE}" ]; then
            UPLOAD_JOB=$(submit_with_retry \
                --account=owner-guest \
                --partition=kingspeak-guest \
                --time=1:00:00 \
                --mem=4G \
                --dependency="afterany:${MERGE_JOB}" \
                --job-name=upload_E${ENERGY_STR}_s${SEED} \
                --output="$HOME/csv_logs/upload_E${ENERGY_STR}_s${SEED}.out" \
                --error="$HOME/csv_logs/upload_E${ENERGY_STR}_s${SEED}.err" \
                --wrap="
source /uufs/chpc.utah.edu/common/home/u1520754/miniconda3/etc/profile.d/conda.sh
conda activate jupyter_env
python3 ${UPLOAD_SCRIPT} \
    --csv ${FINAL_CSV} \
    --base-path ${BASE_PATH} \
    --drive-dest ${UPLOAD_TRIGGERED_TO_DRIVE} \
    --threshold ${UPLOAD_TRIGGERED_THRESHOLD}
")
            if [ -n "${UPLOAD_JOB}" ]; then
                echo "  Upload job: ${UPLOAD_JOB} -> ${UPLOAD_TRIGGERED_TO_DRIVE}"
            else
                echo "  [E=${ENERGY_STR} s${SEED}] ERROR: Failed to submit upload job"
            fi
        fi
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
