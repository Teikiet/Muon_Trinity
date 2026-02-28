#!/bin/bash
# =============================================================================
# run_corsika8_trinity_chain.slurm
#
# If TEL_X and TEL_Z are BOTH non-zero, skip CORSIKA 8 and reuse the
# cherenkov_hits_base.dat from the (TEL_X=0, TEL_Z=0) run, applying only
# the particle_location_correction.py with the x/z offset.
#
# If TEL_X=0 or TEL_Z=0 (at least one is zero), run the full chain:
#   CORSIKA 8 → save cherenkov_hits_base.dat → correction → container
# =============================================================================

#SBATCH --account=owner-guest
#SBATCH --partition=kingspeak-guest
#SBATCH --time=5:00:00
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

# The "base" geometry tag: same params but with x=0, z=0
BASE_GEOM_TAG="Tilt_pdg${PDG}_E${ENERGY}_zen${ZENITH}_az${AZIMUTH}_h${INJ_HEIGHT}_x0_y0_z0_r${TEL_RADIUS}_s${SEED}"
BASE_SIM_DIR="${OUTPUT_BASE_DIR}/${BASE_GEOM_TAG}"

# The base .dat file lives inside the base run's CORSIKA8 (or CHASM) directory
# After flattening/renaming it ends up here:
BASE_DAT_CORSIKA8="${BASE_SIM_DIR}/CORSIKA8/cherenkov_hits_base.dat"
# During the run (before rename) it might still be here:
BASE_DAT_CHASM="${BASE_SIM_DIR}/CHASM/${BASE_GEOM_TAG}/cherenkov_hits_base.dat"

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
# SECTION 1a: DETERMINE RUN MODE
#
# IS_BASE_RUN=true  → TEL_X == 0 AND TEL_Z == 0 → full CORSIKA 8 run
# IS_BASE_RUN=false → TEL_X != 0 OR TEL_Z != 0  → reuse base, skip CORSIKA 8
# =============================================================================

# For the base run (x=0, z=0), check if output already exists
IS_BASE_RUN=false
if [ "$(echo "${TEL_X} == 0 && ${TEL_Z} == 0" | bc -l)" -eq 1 ]; then
    IS_BASE_RUN=true
fi

echo ""
echo "============================================================"
echo "CORSIKA 8 + Trinity Chain"
echo "SLURM Job ID       : ${SLURM_JOB_ID}"
echo "Geometry tag       : ${GEOM_TAG}"
echo "Base geometry tag  : ${BASE_GEOM_TAG}"
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
echo "Run mode           : $(${IS_BASE_RUN} && echo 'FULL RUN (base x=0, z=0)' || echo 'REUSE BASE (offset x=${TEL_X}, z=${TEL_Z})')" 
echo "Is base run (x=0,z=0): ${IS_BASE_RUN}"
echo "Started at         : $(date)"
echo "============================================================"

# =============================================================================
# SECTION 1b: CHECK IF BASE RUN ALREADY EXISTS (for x=0, z=0 only)
# =============================================================================

if ${IS_BASE_RUN}; then
    # Check if the final output already exists
    FINAL_CARE_ROOT="${SIM_DIR}/CARE/cherenkov_hits.root"
    if [ -f "${FINAL_CARE_ROOT}" ] && [ -f "${BASE_DAT_CORSIKA8}" ]; then
        echo ""
        echo "--- Base run output already exists ---"
        echo "  CARE output : ${FINAL_CARE_ROOT}"
        echo "  Base dat    : ${BASE_DAT_CORSIKA8}"
        echo "  Skipping entire run."
        echo "============================================================"
        exit 0
    fi
fi
# =============================================================================
# SECTION 1b2: CHECK IF REUSE RUN ALREADY EXISTS
# =============================================================================

if ! ${IS_BASE_RUN}; then
    FINAL_CARE_ROOT="${SIM_DIR}/CARE/cherenkov_hits.root"
    if [ -f "${FINAL_CARE_ROOT}" ]; then
        echo ""
        echo "--- Reuse run output already exists ---"
        echo "  CARE output : ${FINAL_CARE_ROOT}"
        echo "  Skipping entire run."
        echo "============================================================"
        exit 0
    fi
fi
# =============================================================================
# SECTION 1c: CLEAN UP PREVIOUS RUN OUTPUT (if any)
# =============================================================================

echo ""
echo "--- SECTION 1c: Cleaning up previous run output ---"
echo ""

if [ -d "${SIM_DIR}" ]; then
    echo "Previous SIM_DIR exists: ${SIM_DIR}"

    for OLD_DIR in \
        "${SIM_DIR}/CHASM" \
        "${SIM_DIR}/CORSIKA8" \
        "${SIM_DIR}/CARE" \
        "${SIM_DIR}/GROPT" \
        "${SIM_DIR}/CIO" \
        "${SIM_DIR}/corsika8_output"; do
        if [ -d "${OLD_DIR}" ]; then
            rm -rf "${OLD_DIR}"
            echo "  Removed: ${OLD_DIR}"
        fi
    done

    for OLD_FILE in \
        "${SIM_DIR}/corsika8_tables.dat"; do
        if [ -f "${OLD_FILE}" ]; then
            rm -f "${OLD_FILE}"
            echo "  Removed: ${OLD_FILE}"
        fi
    done

    echo "Cleanup complete."
else
    echo "No previous output found. Clean start."
fi

# =============================================================================
# SECTION 2: SETUP DIRECTORIES AND LOG
# =============================================================================

mkdir -p "${CHASM_TILT_DIR}"
mkdir -p "$(dirname "$LOG_FILE")"

touch "$LOG_FILE"
chmod 666 "$LOG_FILE"

# =============================================================================
# SECTION 3: CORSIKA 8 SIMULATION (or reuse base)
# =============================================================================

if ! ${IS_BASE_RUN}; then
    # -----------------------------------------------------------------
    # REUSE MODE: TEL_X != 0 OR TEL_Z != 0
    # Find and copy the base .dat file instead of running CORSIKA 8
    # -----------------------------------------------------------------
    echo ""
    echo "--- SECTION 3: REUSE MODE — Copying cherenkov_hits_base.dat ---"
    echo "    Offset: x=${TEL_X}, z=${TEL_Z}"
    echo ""

    # Look for the base dat file in possible locations
    BASE_DAT_FOUND=""
    for CANDIDATE in \
        "${BASE_DAT_CORSIKA8}" \
        "${BASE_DAT_CHASM}" \
        "${BASE_SIM_DIR}/CHASM/cherenkov_hits_base.dat" \
        "${BASE_SIM_DIR}/CORSIKA8/${BASE_GEOM_TAG}/cherenkov_hits_base.dat"; do
        if [ -f "${CANDIDATE}" ]; then
            BASE_DAT_FOUND="${CANDIDATE}"
            break
        fi
    done

    if [ -z "${BASE_DAT_FOUND}" ]; then
        echo "ERROR: cherenkov_hits_base.dat not found for base geometry."
        echo "  Expected at: ${BASE_DAT_CORSIKA8}"
        echo "  Or at:       ${BASE_DAT_CHASM}"
        echo ""
        echo "  The base run (TEL_X=0, TEL_Z=0) must complete first."
        echo "  Base SIM_DIR: ${BASE_SIM_DIR}"
        exit 1
    fi

    echo "Found base dat file: ${BASE_DAT_FOUND}"
    echo "  Size: $(du -sh "${BASE_DAT_FOUND}" | cut -f1)"

    # Copy base dat to our CHASM tilt directory as the starting point
    cp "${BASE_DAT_FOUND}" "${CHASM_TILT_DIR}/cherenkov_hits.dat"
    echo "Copied to: ${CHASM_TILT_DIR}/cherenkov_hits.dat"

else
    # -----------------------------------------------------------------
    # FULL RUN MODE: TEL_X == 0 AND TEL_Z == 0 (base run)
    # Run CORSIKA 8 from scratch
    # -----------------------------------------------------------------
    echo ""
    echo "--- SECTION 3: FULL RUN — Running CORSIKA 8 (base geometry) ---"
    echo ""

    mkdir -p "${C8_WORK_DIR}"

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
        --telescope-x 0 \
        --telescope-y 0 \
        --telescope-z 0 \
        --telescope-radius "${CHERENKOV_RADIUS}" \
        --seed "${SEED}" \
        --hadronModel "${HADRON_MODEL}" \
        -f "${C8_OUT_PREFIX}" \
        --cherenkov-output-format eventio \
        --atmosphere-absorption-file /uufs/chpc.utah.edu/common/home/u1520754/corsika/modules/data/CHERENKOV/atmosphere/abstable_MODTRAN_new.dat

    C8_EXIT=$?

    if [ ${C8_EXIT} -ne 0 ]; then
        echo "ERROR: CORSIKA 8 failed with exit code ${C8_EXIT}"
        rm -rf "${C8_WORK_DIR}"
        exit ${C8_EXIT}
    fi

    echo ""
    echo "CORSIKA 8 completed successfully."
    ls -lh "${C8_WORK_DIR}/"

    # -----------------------------------------------------------------
    # Move output files
    # -----------------------------------------------------------------
    echo ""
    echo "--- SECTION 3a: Moving CORSIKA 8 output files ---"
    echo ""

    if [ ! -f "${C8_DAT_FILE}" ]; then
        echo "ERROR: Expected CORSIKA 8 output file not found: ${C8_DAT_FILE}"
        ls -lh "${C8_WORK_DIR}/" || true
        rm -rf "${C8_WORK_DIR}"
        exit 1
    fi

    echo "Found: ${C8_DAT_FILE}  ($(du -sh "${C8_DAT_FILE}" | cut -f1))"

    # Save the FULL cherenkov radius output as cherenkov_hits_base.dat
    cp "${C8_DAT_FILE}" "${CHASM_TILT_DIR}/cherenkov_hits_base.dat"
    echo "Saved base dat:  ${CHASM_TILT_DIR}/cherenkov_hits_base.dat"

    # Also place it as cherenkov_hits.dat (will be filtered by correction script)
    mv "${C8_DAT_FILE}" "${CHASM_TILT_DIR}/cherenkov_hits.dat"
    echo "Moved dat to:    ${CHASM_TILT_DIR}/cherenkov_hits.dat"

    # Move the parquet/yaml output subdir to corsika8_output/
    C8_OUTPUT_SUBDIR="${C8_WORK_DIR}/${C8_OUT_PREFIX}"
    if [ -d "${C8_OUTPUT_SUBDIR}" ]; then
        if [ -d "${SIM_DIR}/corsika8_output" ]; then
            rm -rf "${SIM_DIR}/corsika8_output"
        fi
        mv "${C8_OUTPUT_SUBDIR}" "${SIM_DIR}/corsika8_output"
        echo "CORSIKA 8 output dir moved to: ${SIM_DIR}/corsika8_output"
    fi

    # Move tables.dat if present
    if [ -f "${C8_WORK_DIR}/tables.dat" ]; then
        mv "${C8_WORK_DIR}/tables.dat" "${SIM_DIR}/corsika8_tables.dat"
        echo "Moved tables.dat to: ${SIM_DIR}/corsika8_tables.dat"
    fi

    # Clean up work dir
    cd "${SIM_DIR}"
    rm -rf "${C8_WORK_DIR}"
    echo "Removed work dir: ${C8_WORK_DIR}"
fi


echo ""
echo "CHASM tilt directory now contains:"
ls -lh "${CHASM_TILT_DIR}/"

# =============================================================================
# SECTION 4: PARTICLE LOCATION CORRECTION (muon magnetic field deflection)
#
# For FULL RUN mode (at least one of x,y is 0):
#   - Finds the Cherenkov hotspot peak
#   - Filters by TEL_RADIUS around (peak + offset)
#   - Recenters to origin
#
# For REUSE mode (both x,z non-zero):
#   - Applies TEL_X, TEL_Z offset to the base dat
#   - Filters by TEL_RADIUS around offset position
#   - Recenters to origin
# =============================================================================

echo ""
echo "--- SECTION 4: Particle location correction ---"
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
echo "Telescope offset:  x=${TEL_X}, y=${TEL_Z}"
echo "Recenter:          yes"

# Back up the original file (only if not already backed up from base)
if [ ! -f "${CHERENKOV_DAT}.orig" ]; then
    cp "${CHERENKOV_DAT}" "${CHERENKOV_DAT}.orig"
    echo "Backed up original to: ${CHERENKOV_DAT}.orig"
fi

# Run the correction with the telescope offset
python3 "${CORRECTION_SCRIPT}" \
    --cherenkov-eventio "${CHERENKOV_DAT}" \
    --telescope-radius "${TEL_RADIUS}" \
    --telescope-x "${TEL_X}" \
    --telescope-y "${TEL_Z}" \
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

# Move cherenkov_hits_base.dat OUT of the CHASM tilt dir so GROPT doesn't read it
BASE_DAT_IN_CHASM="${CHASM_TILT_DIR}/cherenkov_hits_base.dat"
BASE_DAT_TEMP="${SIM_DIR}/cherenkov_hits_base.dat.tmp"

if [ -f "${BASE_DAT_IN_CHASM}" ]; then
    mv "${BASE_DAT_IN_CHASM}" "${BASE_DAT_TEMP}"
    echo "Temporarily moved cherenkov_hits_base.dat out of CHASM tilt dir"
    echo "  From: ${BASE_DAT_IN_CHASM}"
    echo "  To:   ${BASE_DAT_TEMP}"
fi

# Also move the .orig backup out if present
ORIG_DAT_IN_CHASM="${CHASM_TILT_DIR}/cherenkov_hits.dat.orig"
ORIG_DAT_TEMP="${SIM_DIR}/cherenkov_hits.dat.orig.tmp"

if [ -f "${ORIG_DAT_IN_CHASM}" ]; then
    mv "${ORIG_DAT_IN_CHASM}" "${ORIG_DAT_TEMP}"
    echo "Temporarily moved cherenkov_hits.dat.orig out of CHASM tilt dir"
fi

chmod -R 777 "${SIM_DIR}"
chmod 666 "${LOG_FILE}"

# =============================================================================
# SECTION 6: CONTAINER — CIO READER → GROPTOCS → CARE
# =============================================================================

echo ""
echo "--- SECTION 6: Running container chain (CIO → GrOptics → CARE) ---"
echo ""

cd "${SIM_DIR}"

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

# Move cherenkov_hits_base.dat back into the CHASM tilt dir
if [ -f "${BASE_DAT_TEMP}" ]; then
    mv "${BASE_DAT_TEMP}" "${BASE_DAT_IN_CHASM}"
    echo "Restored cherenkov_hits_base.dat to: ${BASE_DAT_IN_CHASM}"
fi

if [ -f "${ORIG_DAT_TEMP}" ]; then
    mv "${ORIG_DAT_TEMP}" "${ORIG_DAT_IN_CHASM}"
    echo "Restored cherenkov_hits.dat.orig to: ${ORIG_DAT_IN_CHASM}"
fi

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

cd "${SIM_DIR}"

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
