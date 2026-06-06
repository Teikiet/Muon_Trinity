#!/bin/bash
#SBATCH --job-name=submitter
#SBATCH --time=70:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=kingspeak-guest
#SBATCH --account=owner-guest
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1520754/logs/submitter_%j.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/u1520754/logs/submitter_%j.err

set -eo pipefail

SIM_INPUT_SCRIPT="$HOME/Muon_Trinity/cluster_corsika8/Trinity_sim_imput.py"
MAX_TOTAL_JOBS=500
JOB_POLL_INTERVAL=120

wait_for_jobs_with_prefix() {
    local prefix="$1"
    local poll_interval="$2"
    while true; do
        REMAINING=$(squeue -u "$USER" -h -o "%j" | awk -v p="$prefix" 'index($0,p)==1 {c++} END{print c+0}')
        if [ "${REMAINING}" -eq 0 ]; then
            echo "All jobs with prefix '${prefix}' finished."
            break
        fi
        echo "$(date): ${REMAINING} jobs with prefix '${prefix}' still queued/running..."
        sleep "${poll_interval}"
    done
}

wait_for_total_jobs_below() {
    local limit="$1"
    local poll_interval="$2"
    while true; do
        NJOBS=$(squeue -u "$USER" -h 2>/dev/null | wc -l)
        if [ "${NJOBS}" -lt "${limit}" ]; then
            return
        fi
        echo "$(date): Queue at ${NJOBS}/${limit} jobs. Waiting ${poll_interval}s..."
        sleep "${poll_interval}"
    done
}

source /uufs/chpc.utah.edu/common/home/u1520754/miniconda3/etc/profile.d/conda.sh
conda activate jupyter_env

SIM_INPUT_JSON=$(mktemp /tmp/trinity_sim_input_XXXXXX.json)
trap 'rm -f "${SIM_INPUT_JSON}"' EXIT
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

echo "Found ${#ENERGY_STRS[@]} energies"
echo "Energies: ${ENERGY_STRS[*]}"
echo "Seeds: ${SEEDS[*]}"

CSV_PIDS=()

for ENERGY_STR in "${ENERGY_STRS[@]}"; do
    for SEED in "${SEEDS[@]}"; do
        echo "============================================================"
        echo "ENERGY = ${ENERGY_STR}, SEED = ${SEED}"
        echo "============================================================"

        wait_for_total_jobs_below "${MAX_TOTAL_JOBS}" "${JOB_POLL_INTERVAL}"

        # Phase 1: submit simulation jobs for this energy/seed
        bash $HOME/Muon_Trinity/cluster_corsika8/base_muon_event_generation/submit_multiple_continuous_tree.sh \
            --update_time_stamp\
            --delete_log_on_success\
            --only-energy "${ENERGY_STR}"\
            --only-seed "${SEED}"
            #--rerun-event all 
        # Phase 2+3: wait for this energy/seed to finish, then run CSV in background
        (
            SIM_JOB_PREFIX="corsika8_trinity_pid${PDG}_E${ENERGY_STR}_s${SEED}"
            echo "Waiting for jobs with name prefix '${SIM_JOB_PREFIX}' to finish..."
            wait_for_jobs_with_prefix "${SIM_JOB_PREFIX}" 120

            wait_for_total_jobs_below "${MAX_TOTAL_JOBS}" "${JOB_POLL_INTERVAL}"

            echo "Starting CSV summary collection for E=${ENERGY_STR} s=${SEED}..."
            bash $HOME/Muon_Trinity/cluster_corsika8/save_data2csv/submit_save_care2csv_tree.sh \
                --only-energy "${ENERGY_STR}" \
                --only-seed "${SEED}" \
                --wait\
                --Max_PE_cut 10 
            echo "Completed CSV summary for E=${ENERGY_STR} s=${SEED}."
        ) &
        CSV_PIDS+=("$!")
    done
done

if [ ${#CSV_PIDS[@]} -gt 0 ]; then
    echo "Waiting for ${#CSV_PIDS[@]} CSV tasks to finish..."
    CSV_FAILED=0
    for pid in "${CSV_PIDS[@]}"; do
        if ! wait "${pid}"; then
            CSV_FAILED=1
        fi
    done
    if [ "${CSV_FAILED}" -ne 0 ]; then
        echo "ERROR: One or more CSV tasks failed."
        exit 1
    fi
fi

echo "Done."
