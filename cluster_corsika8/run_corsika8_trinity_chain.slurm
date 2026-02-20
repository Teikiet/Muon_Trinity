#!/bin/bash
# =============================================================================
# run_corsika8_trinity_chain.slurm
# =============================================================================

#SBATCH --account=owner-guest
#SBATCH --partition=kingspeak-guest
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --job-name=corsika8_trinity
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1520754/log/corsika8_trinity_%j.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/u1520754/log/corsika8_trinity_%j.err

set -euo pipefail

# =============================================================================
# SECTION 0: ARGUMENT PARSING
# =============================================================================

if [ $# -lt 14 ]; then
    echo "ERROR: Not enough arguments."
    echo "Usage: sbatch $0 <output_base_dir> <log_file> <pdg> <energy_eV> \
<zenith_deg> <azimuth_deg> <injection_height_m> <obs_level_m> \
<telescope_x> <telescope_y> <telescope_z> <cherenkov_radius_m> \
<telescope_radius_m> <seed> [hadron_model]"
    exit 1
fi

OUTPUT_BASE_DIR="$1"
LOG_FILE="$2"
PDG="$3"
ENERGY="$4"
ZENITH="$5"
AZIMUTH="$6"
INJ_HEIGHT="$7"
OBS_LEVEL="$8"
TEL_X="$9"
TEL_Y="${10}"
TEL_Z="${11}"
CHERENKOV_RADIUS="${12}"   # large sampling area for CORSIKA 8 photon collection
TEL_RADIUS="${13}"         # actual telescope aperture radius for filtering
SEED="${14}"
HADRON_MODEL="${15:-SIBYLL-2.3d}"

# =============================================================================
# SECTION 1: DERIVED PATHS AND NAMES
# =============================================================================

GEOM_TAG="Tilt_pdg${PDG}_E${ENERGY}_zen${ZENITH}_az${AZIMUTH}_h${INJ_HEIGHT}_x${TEL_X}_y${TEL_Y}_z${TEL_Z}_r${TEL_RADIUS}_s${SEED}"

SIM_DIR="${OUTPUT_BASE_DIR}/${GEOM_TAG}"

# During the run, keep the directory named CHASM so the container finds it.
# We rename it to CORSIKA8 after everything is done.
CHASM_TILT_DIR="${SIM_DIR}/CHASM/${GEOM_TAG}"

C8_WORK_DIR="${OUTPUT_BASE_DIR}/corsika8_work_${SLURM_JOB_ID}"
C8_OUT_PREFIX="cherenkov_hits"
C8_DAT_FILE="${C8_WORK_DIR}/${C8_OUT_PREFIX}.dat"

# Path to the particle location correction script (next to this slurm script)
SCRIPT_DIR="$HOME/Muon_Trinity/cluster_corsika8"
CORRECTION_SCRIPT="${SCRIPT_DIR}/particle_location_correction.py"

SANDBOX_DIR="/uufs/chpc.utah.edu/common/home/u1520754/Detector_Sims/Trinity_fullchain/root"
CHAIN_SCRIPT="/root/TrinitySims/BatchSubmissionScripts/Hive_SequentialChainJobPerShower/runTrinityShowersMCSequentialChain_2.sh"

# =============================================================================
# SECTION 2: SETUP DIRECTORIES AND LOG
# =============================================================================

mkdir -p "${CHASM_TILT_DIR}"
mkdir -p "${C8_WORK_DIR}"
mkdir -p "$(dirname "$LOG_FILE")"

touch "$LOG_FILE"
chmod 666 "$LOG_FILE"

echo "============================================================"
echo "CORSIKA 8 + Trinity Chain"
echo "SLURM Job ID       : ${SLURM_JOB_ID}"
echo "Geometry tag       : ${GEOM_TAG}"
echo "SIM_DIR            : ${SIM_DIR}"
echo "CHASM tilt dir     : ${CHASM_TILT_DIR}"
echo "C8 work dir        : ${C8_WORK_DIR}"
echo "Log file           : ${LOG_FILE}"
echo "PDG                : ${PDG}"
echo "Energy (eV)        : ${ENERGY}"
echo "Zenith (deg)       : ${ZENITH}"
echo "Azimuth (deg)      : ${AZIMUTH}"
echo "Inj. height (m)    : ${INJ_HEIGHT}"
echo "Obs. level (m)     : ${OBS_LEVEL}"
echo "Tel. pos (m)       : (${TEL_X}, ${TEL_Y}, ${TEL_Z})"
echo "Cherenkov radius (m): ${CHERENKOV_RADIUS}"
echo "Telescope radius (m): ${TEL_RADIUS}"
echo "Seed               : ${SEED}"
echo "Hadron model       : ${HADRON_MODEL}"
echo "Correction script  : ${CORRECTION_SCRIPT}"
echo "Started at         : $(date)"
echo "============================================================"

# =============================================================================
# SECTION 3: CORSIKA 8 SIMULATION
# =============================================================================

echo ""
echo "--- SECTION 3: Running CORSIKA 8 ---"
echo ""

source ~/miniconda3/etc/profile.d/conda.sh
conda activate corsika

export PATH=$HOME/corsika-install/bin:$PATH

if ! command -v muon_cherenkov &>/dev/null; then
    echo "ERROR: muon_cherenkov not found in PATH after activating conda env."
    echo "PATH=${PATH}"
    exit 1
fi

cd "${C8_WORK_DIR}"

muon_cherenkov \
    --pdg "${PDG}" \
    -E "${ENERGY}" \
    --injection-height "${INJ_HEIGHT}" \
    --emcut 0.02 --mucut 0.43 --hadcut 0.5 \
    --observation-level "${OBS_LEVEL}" \
    --zenith "${ZENITH}" \
    --azimuth "${AZIMUTH}" \
    --cherenkov-projection telescope \
    --telescope-azimuth 270 \
    --telescope-zenith 90 \
    --telescope-x "${TEL_X}" \
    --telescope-y "${TEL_Y}" \
    --telescope-z "${TEL_Z}" \
    --telescope-radius "${CHERENKOV_RADIUS}" \
    --seed "${SEED}" \
    --hadronModel "${HADRON_MODEL}" \
    -f "${C8_OUT_PREFIX}" \
    --cherenkov-output-format eventio \
    --atmosphere-absorption-file /uufs/chpc.utah.edu/common/home/u1520754/corsika/modules/data/CHERENKOV/atmosphere/abstable_MODTRAN_new.dat

C8_EXIT=$?

if [ ${C8_EXIT} -ne 0 ]; then
    echo "ERROR: CORSIKA 8 failed with exit code ${C8_EXIT}"
    exit ${C8_EXIT}
fi

echo ""
echo "CORSIKA 8 completed successfully."
ls -lh "${C8_WORK_DIR}/"

# =============================================================================
# SECTION 4a: MOVE .dat FILE INTO CHASM TILT DIRECTORY
# =============================================================================

echo ""
echo "--- SECTION 4a: Moving .dat file into CHASM tilt directory ---"
echo ""

if [ ! -f "${C8_DAT_FILE}" ]; then
    echo "ERROR: Expected CORSIKA 8 output file not found: ${C8_DAT_FILE}"
    ls -lh "${C8_WORK_DIR}/" || true
    exit 1
fi

echo "Found: ${C8_DAT_FILE}  ($(du -sh "${C8_DAT_FILE}" | cut -f1))"

mv "${C8_DAT_FILE}" "${CHASM_TILT_DIR}/cherenkov_hits.dat"
echo "Moved .dat to: ${CHASM_TILT_DIR}/cherenkov_hits.dat"

# Move the parquet/yaml output subdir to corsika8_output/
C8_OUTPUT_SUBDIR="${C8_WORK_DIR}/${C8_OUT_PREFIX}"
if [ -d "${C8_OUTPUT_SUBDIR}" ]; then
    mv "${C8_OUTPUT_SUBDIR}" "${SIM_DIR}/corsika8_output"
    echo "CORSIKA 8 output dir moved to: ${SIM_DIR}/corsika8_output"
fi

# Move tables.dat if present
if [ -f "${C8_WORK_DIR}/tables.dat" ]; then
    mv "${C8_WORK_DIR}/tables.dat" "${SIM_DIR}/corsika8_tables.dat"
    echo "Moved tables.dat to: ${SIM_DIR}/corsika8_tables.dat"
fi

# Clean up the now-empty work dir
rm -rf "${C8_WORK_DIR}"
echo "Removed work dir: ${C8_WORK_DIR}"

echo ""
echo "CHASM tilt directory now contains:"
ls -lh "${CHASM_TILT_DIR}/"

# =============================================================================
# SECTION 4b: PARTICLE LOCATION CORRECTION (muon magnetic field deflection)
#
# Uses two-pass histogram peak-finding to locate the actual Cherenkov photon
# hotspot (which may be offset from the telescope center due to magnetic
# deflection of the muon), then filters by the *telescope* radius around
# that peak and recenters photon coordinates to the origin.
#
# CHERENKOV_RADIUS was used in CORSIKA 8 to collect photons over a wide area.
# TEL_RADIUS is the actual telescope aperture used here for the final filter.
#
# This overwrites cherenkov_hits.dat in place.
# =============================================================================

echo ""
echo "--- SECTION 4b: Particle location correction (magnetic deflection) ---"
echo ""

CHERENKOV_DAT="${CHASM_TILT_DIR}/cherenkov_hits.dat"

if [ ! -f "${CORRECTION_SCRIPT}" ]; then
    echo "ERROR: Correction script not found: ${CORRECTION_SCRIPT}"
    exit 1
fi

if [ ! -f "${CHERENKOV_DAT}" ]; then
    echo "ERROR: cherenkov_hits.dat not found at: ${CHERENKOV_DAT}"
    exit 1
fi

echo "Input file:        ${CHERENKOV_DAT}  ($(du -sh "${CHERENKOV_DAT}" | cut -f1))"
echo "Cherenkov radius:  ${CHERENKOV_RADIUS} m  (CORSIKA 8 sampling area)"
echo "Telescope radius:  ${TEL_RADIUS} m  (final aperture filter)"
echo "Recenter:          yes"

# Back up the original file
cp "${CHERENKOV_DAT}" "${CHERENKOV_DAT}.orig"
echo "Backed up original to: ${CHERENKOV_DAT}.orig"

# Run the correction: filter around the peak and recenter
python3 "${CORRECTION_SCRIPT}" \
    --cherenkov-eventio "${CHERENKOV_DAT}" \
    --telescope-radius "${TEL_RADIUS}" \
    --output "${CHERENKOV_DAT}.corrected" \
    --recenter \
    --nbins 100

CORR_EXIT=$?

if [ ${CORR_EXIT} -ne 0 ]; then
    echo "ERROR: Particle location correction failed with exit code ${CORR_EXIT}"
    echo "Restoring original file."
    mv "${CHERENKOV_DAT}.orig" "${CHERENKOV_DAT}"
    exit ${CORR_EXIT}
fi

# Replace the original with the corrected version
mv "${CHERENKOV_DAT}.corrected" "${CHERENKOV_DAT}"
echo ""
echo "Corrected file written: ${CHERENKOV_DAT}  ($(du -sh "${CHERENKOV_DAT}" | cut -f1))"
echo "Original backup kept:   ${CHERENKOV_DAT}.orig"

# =============================================================================
# SECTION 5: PERMISSIONS BEFORE CONTAINER
# =============================================================================

echo ""
echo "--- SECTION 5: Setting permissions ---"
chmod -R 777 "${SIM_DIR}"
chmod 666 "${LOG_FILE}"

# =============================================================================
# SECTION 6: CONTAINER — CIO READER → GROPTOCS → CARE
# =============================================================================

echo ""
echo "--- SECTION 6: Running container chain (CIO → GrOptics → CARE) ---"
echo ""

module load apptainer

LOG_FILE_DIR="$(dirname "${LOG_FILE}")"
LOG_FILE_NAME="$(basename "${LOG_FILE}")"

apptainer exec --pid --cleanenv \
    --bind "${SIM_DIR}:/in" \
    --bind "${LOG_FILE_DIR}:/output_mnt" \
    "${SANDBOX_DIR}" \
    bash "${CHAIN_SCRIPT}" \
        "/in" \
        "/output_mnt/${LOG_FILE_NAME}" \
        "${SEED}"

CONTAINER_EXIT=$?

echo ""
echo "Container chain exit code: ${CONTAINER_EXIT}"

# =============================================================================
# SECTION 7: POST-RUN CLEANUP — FLATTEN AND RENAME
# =============================================================================

echo ""
echo "--- SECTION 7: Flattening nested Tilt_ subdirectories ---"
echo ""

flatten_tilt_subdir() {
    local PARENT_DIR="$1"
    local LABEL="$2"

    if [ ! -d "${PARENT_DIR}" ]; then
        echo "  [${LABEL}] Directory not found, skipping: ${PARENT_DIR}"
        return
    fi

    local FOUND=0
    for SUBDIR in "${PARENT_DIR}"/Tilt_*/; do
        [ -d "${SUBDIR}" ] || continue
        FOUND=1

        echo "  [${LABEL}] Flattening: ${SUBDIR} → ${PARENT_DIR}/"

        for ITEM in "${SUBDIR}"*; do
            [ -e "${ITEM}" ] || continue
            ITEM_NAME="$(basename "${ITEM}")"
            DEST="${PARENT_DIR}/${ITEM_NAME}"
            if [ -e "${DEST}" ]; then
                echo "    WARNING: ${DEST} already exists, skipping ${ITEM_NAME}"
            else
                mv "${ITEM}" "${DEST}"
                echo "    Moved: ${ITEM_NAME}"
            fi
        done

        rmdir "${SUBDIR}" 2>/dev/null \
            && echo "  [${LABEL}] Removed empty subdir: $(basename "${SUBDIR}")" \
            || echo "  [${LABEL}] WARNING: Could not remove ${SUBDIR} (not empty?)"
    done

    if [ ${FOUND} -eq 0 ]; then
        echo "  [${LABEL}] No Tilt_* subdirectory found, nothing to flatten."
    fi
}

flatten_tilt_subdir "${SIM_DIR}/CHASM" "CHASM"
flatten_tilt_subdir "${SIM_DIR}/CARE"  "CARE"
flatten_tilt_subdir "${SIM_DIR}/GROPT" "GROPT"
flatten_tilt_subdir "${SIM_DIR}/CIO"   "CIO"

echo ""
echo "--- Renaming CHASM/ → CORSIKA8/ ---"
if [ -d "${SIM_DIR}/CHASM" ]; then
    mv "${SIM_DIR}/CHASM" "${SIM_DIR}/CORSIKA8"
    echo "Renamed: ${SIM_DIR}/CHASM → ${SIM_DIR}/CORSIKA8"
else
    echo "WARNING: ${SIM_DIR}/CHASM not found, skipping rename."
fi

# =============================================================================
# SECTION 8: FINAL STATUS
# =============================================================================

echo ""
echo "============================================================"
if [ ${CONTAINER_EXIT} -eq 0 ]; then
    echo "SUCCESS: Full chain completed."
    echo ""
    echo "Final output structure:"
    find "${SIM_DIR}" -maxdepth 4 -type f | sort
else
    echo "ERROR: Container chain failed with exit code ${CONTAINER_EXIT}."
fi
echo "Finished at: $(date)"
echo "============================================================"

exit ${CONTAINER_EXIT}
