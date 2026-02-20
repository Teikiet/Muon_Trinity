#!/bin/bash
# =============================================================================
# submit_multiple.sh
#
# Submit multiple run_corsika8_trinity_chain.slurm jobs sweeping over
# azimuth, zenith, telescope x/z, injection height, and seed.
# All outputs go into a single parent directory:
#   <BASE>/<OUTPUT_DIR_NAME>/
#
# Usage:
#   bash submit_multiple.sh
# =============================================================================

# =============================================================================
# FIXED PARAMETERS (same for all jobs)
# =============================================================================
PDG=13                   # muon-
ENERGY=1e4               # primary energy in eV
TEL_Y=0                  # telescope y position (fixed)
CHERENKOV_RADIUS=150     # large sampling area for CORSIKA 8 photon collection (m)
TEL_RADIUS=5             # actual telescope aperture radius (m)
OBS_LEVEL=2944           # observation level in meters
HADRON_MODEL="SIBYLL-2.3d"

SLURM_SCRIPT="$HOME/Muon_Trinity/cluster_corsika8/run_corsika8_trinity_chain.slurm"

# Parent output directory (all jobs share this)
BASE_DIR="$HOME/test_sim_chain"
OUTPUT_DIR_NAME="Muon_pid${PDG}_E${ENERGY}_R${TEL_RADIUS}"
OUTPUT_BASE_DIR="${BASE_DIR}/${OUTPUT_DIR_NAME}"

# Log directory (one log per job, named by geometry tag)
LOG_DIR="${OUTPUT_BASE_DIR}/logs"
mkdir -p "${LOG_DIR}"

# =============================================================================
# VARIABLE PARAMETERS
# Define arrays for each swept variable.
# Every combination of (ZENITH, AZIMUTH, TEL_X, TEL_Z, INJ_HEIGHT, SEED)
# will be submitted as a separate job.
# =============================================================================

ZENITHS=(89.9 89 88 87 86)                  # zenith angles in deg 
AZIMUTHS=(266 267 268 269 270 271 272 273 274)                # azimuth angles in degrees 
TEL_XS=(0)                    # telescope x positions in meters
TEL_ZS=(0)                    # telescope z positions in meters
INJ_HEIGHTS=(5000 10000 15000 20000)  # injection heights in meters
SEEDS=(5)                     # random seeds

# =============================================================================
# SUBMISSION LOOP
# =============================================================================

TOTAL=0
SUBMITTED=0

echo "============================================================"
echo "Batch submission for: ${OUTPUT_DIR_NAME}"
echo "Output base dir     : ${OUTPUT_BASE_DIR}"
echo "SLURM script        : ${SLURM_SCRIPT}"
echo "PDG=${PDG}  ENERGY=${ENERGY}  CHERENKOV_RADIUS=${CHERENKOV_RADIUS}  TEL_RADIUS=${TEL_RADIUS}"
echo "============================================================"
echo ""

for ZENITH in "${ZENITHS[@]}"; do
for AZIMUTH in "${AZIMUTHS[@]}"; do
for TEL_X in "${TEL_XS[@]}"; do
for TEL_Z in "${TEL_ZS[@]}"; do
for INJ_HEIGHT in "${INJ_HEIGHTS[@]}"; do
for SEED in "${SEEDS[@]}"; do

    TOTAL=$((TOTAL + 1))

    # Per-job log file inside the shared log directory
    LOG_FILE="${LOG_DIR}/sim_zen${ZENITH}_az${AZIMUTH}_x${TEL_X}_z${TEL_Z}_h${INJ_HEIGHT}_s${SEED}.log"

    echo "Submitting job ${TOTAL}:"
    echo "  zenith=${ZENITH} azimuth=${AZIMUTH} x=${TEL_X} z=${TEL_Z} h=${INJ_HEIGHT} seed=${SEED}"

    sbatch "${SLURM_SCRIPT}" \
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
        "${HADRON_MODEL}"

    if [ $? -eq 0 ]; then
        SUBMITTED=$((SUBMITTED + 1))
    else
        echo "  WARNING: sbatch failed for this combination."
    fi

    echo ""

done
done
done
done
done
done

# =============================================================================
# SUMMARY
# =============================================================================

echo "============================================================"
echo "Submission complete."
echo "  Total combinations : ${TOTAL}"
echo "  Successfully submitted: ${SUBMITTED}"
echo "  Failed submissions    : $((TOTAL - SUBMITTED))"
echo ""
echo "All outputs will be under:"
echo "  ${OUTPUT_BASE_DIR}/"
echo ""
echo "Per-job SLURM logs (stdout/stderr) are in:"
echo "  $HOME/log/corsika8_trinity_<jobid>.out/err"
echo "Per-job chain logs are in:"
echo "  ${LOG_DIR}/"
echo "============================================================"
