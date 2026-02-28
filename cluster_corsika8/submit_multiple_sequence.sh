#!/bin/bash
# =============================================================================
# submit_multiple_sequence.sh
# Make sure all the based run is submitted before running this script
# Groups jobs into sequences of SEQ_LENGTH. Each sequence is ONE SLURM job
# that runs SEQ_LENGTH simulations back-to-back in a loop.
# Multiple sequences run in parallel across different nodes.
#
# Usage:
#   bash submit_multiple_sequence.sh
# =============================================================================

# =============================================================================
# SEQUENCE LENGTH — how many simulations run back-to-back in one SLURM job
# =============================================================================
SEQ_LENGTH=50 

# Time per single simulation (estimate)
TIME_PER_SIM="1:00:00"
# Total time = SEQ_LENGTH * TIME_PER_SIM (with some padding)
TOTAL_TIME="6:00:00"

# =============================================================================
# FIXED PARAMETERS (same for all jobs)
# =============================================================================
PDG=13
ENERGY=1e5
TEL_Y=0
CHERENKOV_RADIUS=15
TEL_RADIUS=5
OBS_LEVEL=2944
HADRON_MODEL="SIBYLL-2.3d"

# The ORIGINAL single-simulation slurm script (used as a worker)
WORKER_SCRIPT="$HOME/Muon_Trinity/cluster_corsika8/run_corsika8_trinity_chain.slurm"

# Wrapper script that loops over a batch file
WRAPPER_SCRIPT="$HOME/Muon_Trinity/cluster_corsika8/run_sequence_wrapper.slurm"

# Parent output directory
BASE_DIR="/scratch/general/vast/u1520754/muon_sim_chain"
OUTPUT_DIR_NAME="Muon_pid${PDG}_E${ENERGY}_R${TEL_RADIUS}"
OUTPUT_BASE_DIR="${BASE_DIR}/${OUTPUT_DIR_NAME}"

LOG_DIR="${OUTPUT_BASE_DIR}/logs"
BATCH_DIR="${OUTPUT_BASE_DIR}/batch_files"
mkdir -p "${LOG_DIR}" "${BATCH_DIR}"

# =============================================================================
# VARIABLE PARAMETERS 
# position: -10 -7 -5 -2 -1 -0.7 -0.5 -0.2 0 0.2 0.5 0.7 1 2 5 7 10
# production height: 5000 10000 15000 20000 30000
# zenith: 89.9 89 88 87
# azimuth: 268 269 270 271 272
# seed: 5
# =============================================================================

ZENITHS=(89.9 89 88 87)
AZIMUTHS=(268 269 270 271 272)
TEL_XS=(-5 -2 -1 -0.7 -0.5 -0.2 0 0.2 0.5 0.7 1 2 5)
TEL_ZS=(-5 -2 -1 -0.7 -0.5 -0.2 0 0.2 0.5 0.7 1 2 5)
INJ_HEIGHTS=(5000 10000 15000 20000 30000)
SEEDS=(2)

# =============================================================================
# BUILD LIST OF ALL PARAMETER COMBINATIONS
# =============================================================================

declare -a ALL_ARGS
TOTAL=0

for ZENITH in "${ZENITHS[@]}"; do
for AZIMUTH in "${AZIMUTHS[@]}"; do
for TEL_X in "${TEL_XS[@]}"; do
for TEL_Z in "${TEL_ZS[@]}"; do
for INJ_HEIGHT in "${INJ_HEIGHTS[@]}"; do
for SEED in "${SEEDS[@]}"; do

    LOG_FILE="${LOG_DIR}/sim_zen${ZENITH}_az${AZIMUTH}_x${TEL_X}_z${TEL_Z}_h${INJ_HEIGHT}_s${SEED}.log"

    # Store all arguments as a single pipe-delimited line
    ALL_ARGS[$TOTAL]="${OUTPUT_BASE_DIR}|${LOG_FILE}|${PDG}|${ENERGY}|${ZENITH}|${AZIMUTH}|${INJ_HEIGHT}|${OBS_LEVEL}|${TEL_X}|${TEL_Y}|${TEL_Z}|${CHERENKOV_RADIUS}|${TEL_RADIUS}|${SEED}|${HADRON_MODEL}"
    TOTAL=$((TOTAL + 1))

done
done
done
done
done
done

NUM_SEQUENCES=$(( (TOTAL + SEQ_LENGTH - 1) / SEQ_LENGTH ))

echo "============================================================"
echo "Batch submission for: ${OUTPUT_DIR_NAME}"
echo "Output base dir     : ${OUTPUT_BASE_DIR}"
echo "PDG=${PDG}  ENERGY=${ENERGY}  CHERENKOV_RADIUS=${CHERENKOV_RADIUS}  TEL_RADIUS=${TEL_RADIUS}"
echo ""
echo "Total simulations    : ${TOTAL}"
echo "Sequence length      : ${SEQ_LENGTH}"
echo "Number of SLURM jobs : ${NUM_SEQUENCES}"
echo "Time per SLURM job   : ${TOTAL_TIME}"
echo "============================================================"
echo ""

# =============================================================================
# WRITE BATCH FILES AND SUBMIT
# =============================================================================

SUBMITTED=0
FAILED=0

for (( seq=0; seq<NUM_SEQUENCES; seq++ )); do

    START_IDX=$(( seq * SEQ_LENGTH ))
    END_IDX=$(( START_IDX + SEQ_LENGTH ))
    if [ ${END_IDX} -gt ${TOTAL} ]; then
        END_IDX=${TOTAL}
    fi
    COUNT=$(( END_IDX - START_IDX ))

    # Write a batch file listing all simulations for this sequence
    BATCH_FILE="${BATCH_DIR}/sequence_$(printf '%04d' ${seq}).txt"
    > "${BATCH_FILE}"
    for (( i=START_IDX; i<END_IDX; i++ )); do
        echo "${ALL_ARGS[$i]}" >> "${BATCH_FILE}"
    done

    SEQ_LOG="${LOG_DIR}/sequence_$(printf '%04d' ${seq}).log"

    echo "Submitting sequence $((seq+1))/${NUM_SEQUENCES} (${COUNT} sims): ${BATCH_FILE}"

    RESULT=$( sbatch \
        --time="${TOTAL_TIME}" \
        --job-name="corsika8_seq$(printf '%04d' ${seq})" \
        --output="/uufs/chpc.utah.edu/common/home/u1520754/log/corsika8_seq$(printf '%04d' ${seq})_%j.out" \
        --error="/uufs/chpc.utah.edu/common/home/u1520754/log/corsika8_seq$(printf '%04d' ${seq})_%j.err" \
        "${WRAPPER_SCRIPT}" \
        "${BATCH_FILE}" \
        "${WORKER_SCRIPT}" 2>&1 )

    if [ $? -eq 0 ]; then
        SUBMITTED=$((SUBMITTED + 1))
        echo "  ${RESULT}"
    else
        echo "  WARNING: sbatch failed: ${RESULT}"
        FAILED=$((FAILED + 1))
    fi
    echo ""

done

# =============================================================================
# SUMMARY
# =============================================================================

echo "============================================================"
echo "Submission complete."
echo "  Total simulations      : ${TOTAL}"
echo "  SLURM jobs submitted   : ${SUBMITTED}"
echo "  Failed submissions     : ${FAILED}"
echo "  Sims per job           : up to ${SEQ_LENGTH}"
echo "  Time per job           : ${TOTAL_TIME}"
echo ""
echo "Batch files in : ${BATCH_DIR}/"
echo "Outputs in     : ${OUTPUT_BASE_DIR}/"
echo "SLURM logs in  : $HOME/log/"
echo "Chain logs in  : ${LOG_DIR}/"
echo "============================================================"
