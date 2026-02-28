#!/bin/bash
# =============================================================================
# submit_multiple_continuous.sh
#
# Submit multiple run_corsika8_trinity_chain.slurm jobs sweeping over
# azimuth, zenith, telescope x/z, injection height, and seed.
#
# Instead of batch dependencies, this script continuously submits jobs
# as long as the queue has room (below MAX_QUEUED). No jobs wait
# unnecessarily — the queue stays maximally full.
#
# Usage:
#   bash submit_multiple_continuous.sh
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

SLURM_SCRIPT="$HOME/Muon_Trinity/cluster_corsika8/run_corsika8_trinity_chain.slurm"

# Parent output directory (all jobs share this)
BASE_DIR="/scratch/general/vast/u1520754/muon_sim_chain"
OUTPUT_DIR_NAME="Muon_pid${PDG}_E${ENERGY}_R${TEL_RADIUS}"
OUTPUT_BASE_DIR="${BASE_DIR}/${OUTPUT_DIR_NAME}"

# Log directory (one log per job, named by geometry tag)
LOG_DIR="${OUTPUT_BASE_DIR}/logs"
mkdir -p "${LOG_DIR}"

# =============================================================================
# QUEUE CONTROL
# =============================================================================
MAX_QUEUED=500    # maximum jobs in queue at any time
POLL_INTERVAL=30  # seconds between queue checks when full

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

ZENITHS=(89.9 89 88 87)
AZIMUTHS=(268 269 270 271 272)
TEL_XS=(0)
TEL_ZS=(0)
INJ_HEIGHTS=(5000 10000 15000 20000 30000)
SEEDS=($(seq 9 100))

# =============================================================================
# SUBMISSION LOOP — continuous fill, no dependencies
# =============================================================================

TOTAL=0
SUBMITTED=0
SKIPPED=0

echo "============================================================"
echo "Continuous submission for: ${OUTPUT_DIR_NAME}"
echo "Output base dir     : ${OUTPUT_BASE_DIR}"
echo "SLURM script        : ${SLURM_SCRIPT}"
echo "PDG=${PDG}  ENERGY=${ENERGY}  CHERENKOV_RADIUS=${CHERENKOV_RADIUS}  TEL_RADIUS=${TEL_RADIUS}"
echo "Max queued jobs     : ${MAX_QUEUED}"
echo "Seeds               : ${SEEDS[0]} to ${SEEDS[-1]} (${#SEEDS[@]} total)"
echo "============================================================"
echo ""

for SEED in "${SEEDS[@]}"; do
for ZENITH in "${ZENITHS[@]}"; do
for AZIMUTH in "${AZIMUTHS[@]}"; do
for TEL_X in "${TEL_XS[@]}"; do
for TEL_Z in "${TEL_ZS[@]}"; do
for INJ_HEIGHT in "${INJ_HEIGHTS[@]}"; do

    TOTAL=$((TOTAL + 1))

    # -------------------------------------------------------------------------
    # SKIP if output already exists (resume-friendly)
    # -------------------------------------------------------------------------
    GEOM_TAG="Tilt_zen${ZENITH}_az${AZIMUTH}_x${TEL_X}_z${TEL_Z}_h${INJ_HEIGHT}_s${SEED}"
    CARE_OUTPUT="${OUTPUT_BASE_DIR}/${GEOM_TAG}/care_output"
    if [ -d "${CARE_OUTPUT}" ] && [ "$(ls -A "${CARE_OUTPUT}" 2>/dev/null)" ]; then
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    LOG_FILE="${LOG_DIR}/sim_zen${ZENITH}_az${AZIMUTH}_x${TEL_X}_z${TEL_Z}_h${INJ_HEIGHT}_s${SEED}.log"

    # Wait until there's room in the queue
    wait_for_queue_room

    # Submit the job (no dependencies)
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

        # Print progress every 50 jobs to reduce spam
        if (( SUBMITTED % 50 == 0 )); then
            NJOBS=$(squeue -u "$USER" -h 2>/dev/null | wc -l)
            echo "  [$(date +%H:%M:%S)] Submitted: ${SUBMITTED} | Queue: ${NJOBS} | Current: s=${SEED} zen=${ZENITH} az=${AZIMUTH} h=${INJ_HEIGHT}"
        fi
    else
        echo "    WARNING: sbatch failed: ${SBATCH_OUTPUT}"

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
            else
                echo "    WARNING: Retry also failed: ${SBATCH_OUTPUT}"
            fi
        fi
    fi

done
done
done
done
done
done

# =============================================================================
# SUMMARY
# =============================================================================

echo ""
echo "============================================================"
echo "Submission complete."
echo "  Total combinations    : ${TOTAL}"
echo "  Already done (skipped): ${SKIPPED}"
echo "  Successfully submitted: ${SUBMITTED}"
echo "  Failed submissions    : $((TOTAL - SUBMITTED - SKIPPED))"
echo ""
echo "All outputs will be under:"
echo "  ${OUTPUT_BASE_DIR}/"
echo ""
echo "Per-job SLURM logs (stdout/stderr) are in:"
echo "  $HOME/log/corsika8_trinity_<jobid>.out/err"
echo "Per-job chain logs are in:"
echo "  ${LOG_DIR}/"
echo "============================================================"
