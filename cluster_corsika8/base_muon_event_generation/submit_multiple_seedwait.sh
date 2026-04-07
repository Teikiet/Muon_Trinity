#!/bin/bash
# =============================================================================
# submit_multiple_seedwait.sh
#
# Submit multiple run_corsika8_trinity_chain.slurm jobs sweeping over
# azimuth, zenith, telescope x/z, injection height, and seed.
# Seeds are submitted in batches (SEEDS_PER_BATCH at a time).
# All jobs within a batch run in parallel. The next batch waits for
# all jobs in the previous batch to finish before starting.
#
# Includes queue throttling to avoid QOSMaxSubmitJobPerUserLimit errors.
#
# Usage:
#   bash submit_multiple_seedwait.sh
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
# QUEUE THROTTLING
# =============================================================================
MAX_QUEUED=900   # stay safely below your cluster's per-user limit

wait_for_queue_room() {
    while true; do
        NJOBS=$(squeue -u "$USER" -h 2>/dev/null | wc -l)
        if [ "${NJOBS}" -lt "${MAX_QUEUED}" ]; then
            break
        fi
        echo "    Queue full (${NJOBS} jobs queued/running). Waiting 600s..."
        sleep 600
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
# SEED BATCHING
# =============================================================================
SEEDS_PER_BATCH=5   # <--- number of seeds to submit simultaneously

# Split SEEDS array into batches
NUM_SEEDS=${#SEEDS[@]}
NUM_BATCHES=$(( (NUM_SEEDS + SEEDS_PER_BATCH - 1) / SEEDS_PER_BATCH ))

echo "============================================================"
echo "Batch submission for: ${OUTPUT_DIR_NAME}"
echo "Output base dir     : ${OUTPUT_BASE_DIR}"
echo "SLURM script        : ${SLURM_SCRIPT}"
echo "PDG=${PDG}  ENERGY=${ENERGY}  CHERENKOV_RADIUS=${CHERENKOV_RADIUS}  TEL_RADIUS=${TEL_RADIUS}"
echo "Total seeds         : ${NUM_SEEDS}"
echo "Seeds per batch     : ${SEEDS_PER_BATCH}"
echo "Number of batches   : ${NUM_BATCHES}"
echo "============================================================"
echo ""

# =============================================================================
# SUBMISSION LOOP — batch-by-batch with dependency
# =============================================================================

TOTAL=0
SUBMITTED=0
PREV_BATCH_JOBIDS=""   # colon-separated job IDs from the previous batch

for (( BATCH=0; BATCH<NUM_BATCHES; BATCH++ )); do

    # Determine which seeds are in this batch
    START_IDX=$(( BATCH * SEEDS_PER_BATCH ))
    BATCH_SEEDS=("${SEEDS[@]:${START_IDX}:${SEEDS_PER_BATCH}}")

    echo "============================================================"
    echo "BATCH $((BATCH + 1)) of ${NUM_BATCHES}  —  Seeds: ${BATCH_SEEDS[*]}"
    if [ -n "${PREV_BATCH_JOBIDS}" ]; then
        echo "  Dependency: waiting for previous batch jobs"
    else
        echo "  No dependency (first batch)"
    fi
    echo "============================================================"

    CURRENT_BATCH_JOBIDS=""   # collect job IDs for this batch

    for SEED in "${BATCH_SEEDS[@]}"; do
    for ZENITH in "${ZENITHS[@]}"; do
    for AZIMUTH in "${AZIMUTHS[@]}"; do
    for TEL_X in "${TEL_XS[@]}"; do
    for TEL_Z in "${TEL_ZS[@]}"; do
    for INJ_HEIGHT in "${INJ_HEIGHTS[@]}"; do

        TOTAL=$((TOTAL + 1))

        LOG_FILE="${LOG_DIR}/sim_zen${ZENITH}_az${AZIMUTH}_x${TEL_X}_z${TEL_Z}_h${INJ_HEIGHT}_s${SEED}.log"

        echo "  Submitting job ${TOTAL}: seed=${SEED} zen=${ZENITH} az=${AZIMUTH} x=${TEL_X} z=${TEL_Z} h=${INJ_HEIGHT}"

        # Wait if queue is near the limit
        wait_for_queue_room

        # Build sbatch command with optional dependency
        SBATCH_CMD="sbatch"
        if [ -n "${PREV_BATCH_JOBIDS}" ]; then
            SBATCH_CMD="${SBATCH_CMD} --dependency=afterany:${PREV_BATCH_JOBIDS}"
        fi

        # Submit and capture job ID
        SBATCH_OUTPUT=$(${SBATCH_CMD} "${SLURM_SCRIPT}" \
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
            echo "    -> Job ID: ${JOBID}"
            SUBMITTED=$((SUBMITTED + 1))

            if [ -n "${CURRENT_BATCH_JOBIDS}" ]; then
                CURRENT_BATCH_JOBIDS="${CURRENT_BATCH_JOBIDS}:${JOBID}"
            else
                CURRENT_BATCH_JOBIDS="${JOBID}"
            fi
        else
            echo "    WARNING: sbatch failed: ${SBATCH_OUTPUT}"
            if echo "${SBATCH_OUTPUT}" | grep -q "QOSMaxSubmitJobPerUserLimit"; then
                echo "    Queue limit hit. Waiting 120s and retrying..."
                sleep 120
                wait_for_queue_room
                SBATCH_OUTPUT=$(${SBATCH_CMD} "${SLURM_SCRIPT}" \
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
                    if [ -n "${CURRENT_BATCH_JOBIDS}" ]; then
                        CURRENT_BATCH_JOBIDS="${CURRENT_BATCH_JOBIDS}:${JOBID}"
                    else
                        CURRENT_BATCH_JOBIDS="${JOBID}"
                    fi
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

    echo ""
    echo "  Batch $((BATCH + 1)): seeds=${BATCH_SEEDS[*]}, jobs submitted for this batch"
    echo ""

    # This batch's jobs become the dependency for the next batch
    PREV_BATCH_JOBIDS="${CURRENT_BATCH_JOBIDS}"

done

# =============================================================================
# SUMMARY
# =============================================================================

echo "============================================================"
echo "Submission complete."
echo "  Total combinations    : ${TOTAL}"
echo "  Successfully submitted: ${SUBMITTED}"
echo "  Failed submissions    : $((TOTAL - SUBMITTED))"
echo "  Seed batches          : ${NUM_BATCHES} (${SEEDS_PER_BATCH} seeds each)"
echo ""
echo "All outputs will be under:"
echo "  ${OUTPUT_BASE_DIR}/"
echo ""
echo "Per-job SLURM logs (stdout/stderr) are in:"
echo "  $HOME/log/corsika8_trinity_<jobid>.out/err"
echo "Per-job chain logs are in:"
echo "  ${LOG_DIR}/"
echo "============================================================"
