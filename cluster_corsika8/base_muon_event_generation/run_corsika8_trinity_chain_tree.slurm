#!/bin/bash
# =============================================================================
# run_corsika8_trinity_chain_tree.slurm
#
# Tree structure:
#   <OUTPUT_BASE_DIR>/pdg${PDG}_E${ENERGY}_r${TEL_RADIUS}_s${SEED}/
#     zen${ZENITH}/
#       az${AZIMUTH}/
#         h${INJ_HEIGHT}/
#           x${TEL_X}_y${TEL_Y}_z${TEL_Z}/
#             CHASM/ (renamed to CORSIKA8/ after run)
#             CARE/
#             GROPT/
#             CIO/
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
#SBATCH --output=/scratch/general/vast/u1520754/log/corsika8_trinity_%j.out
#SBATCH --error=/scratch/general/vast/u1520754/log/corsika8_trinity_%j.err

set -euo pipefail

# =============================================================================
# SECTION 0: ARGUMENT PARSING
# =============================================================================

if [ $# -lt 14 ]; then
    echo "ERROR: Not enough arguments."
    echo "Usage: sbatch $0 <output_base_dir> <log_file> <pdg> <energy_eV> \
<zenith_deg> <azimuth_deg> <injection_height_m> <obs_level_m> \
<telescope_x> <telescope_y> <telescope_z> <cherenkov_radius_m> \
<telescope_radius_m> <seed> [hadron_model] [rerun_event] [--update_time_stamp] [--reduced_based_run_radius[=R]]"
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

# Optional arguments (backward compatible):
#   positional: [hadron_model] [rerun_event]
#   flag: --update_time_stamp (or --update_time_stamp=true/false)
#   flag: --reduced_based_run_radius (or --reduced_based_run_radius=R)
HADRON_MODEL="SIBYLL-2.3d"
RERUN_EVENT_RAW="false"
UPDATE_TIME_STAMP_RAW="false"
REDUCED_BASED_RUN_RADIUS_ENABLED_RAW="false"
REDUCED_BASED_RUN_RADIUS_RAW="15"

POSITIONAL_OPTIONALS=()
if [ $# -ge 15 ]; then
    for ARG in "${@:15}"; do
        case "${ARG}" in
            --update_time_stamp)
                UPDATE_TIME_STAMP_RAW="true"
                ;;
            --update_time_stamp=*)
                UPDATE_TIME_STAMP_RAW="${ARG#*=}"
                ;;
            --hadron_model=*)
                HADRON_MODEL="${ARG#*=}"
                ;;
            --rerun_event=*)
                RERUN_EVENT_RAW="${ARG#*=}"
                ;;
            --reduced_based_run_radius)
                REDUCED_BASED_RUN_RADIUS_ENABLED_RAW="true"
                REDUCED_BASED_RUN_RADIUS_RAW="15"
                ;;
            --reduced_based_run_radius=*)
                REDUCED_BASED_RUN_RADIUS_ENABLED_RAW="true"
                REDUCED_BASED_RUN_RADIUS_RAW="${ARG#*=}"
                ;;
            *)
                POSITIONAL_OPTIONALS+=("${ARG}")
                ;;
        esac
    done
fi

if [ ${#POSITIONAL_OPTIONALS[@]} -ge 1 ]; then
    HADRON_MODEL="${POSITIONAL_OPTIONALS[0]}"
fi
if [ ${#POSITIONAL_OPTIONALS[@]} -ge 2 ]; then
    RERUN_EVENT_RAW="${POSITIONAL_OPTIONALS[1]}"
fi
if [ ${#POSITIONAL_OPTIONALS[@]} -gt 2 ]; then
    echo "ERROR: Too many positional optional arguments."
    echo "Expected at most: [hadron_model] [rerun_event]"
    exit 1
fi

case "${RERUN_EVENT_RAW,,}" in
    true|1|yes|y)
        RERUN_EVENT=true
        ;;
    false|0|no|n|"")
        RERUN_EVENT=false
        ;;
    *)
        echo "ERROR: Invalid rerun_event value '${RERUN_EVENT_RAW}'. Use true/false."
        exit 1
        ;;
esac

case "${UPDATE_TIME_STAMP_RAW,,}" in
    true|1|yes|y)
        UPDATE_TIME_STAMP=true
        ;;
    false|0|no|n|"")
        UPDATE_TIME_STAMP=false
        ;;
    *)
        echo "ERROR: Invalid --update_time_stamp value '${UPDATE_TIME_STAMP_RAW}'. Use true/false."
        exit 1
        ;;
esac

case "${REDUCED_BASED_RUN_RADIUS_ENABLED_RAW,,}" in
    true|1|yes|y)
        REDUCED_BASED_RUN_RADIUS_ENABLED=true
        ;;
    false|0|no|n|"")
        REDUCED_BASED_RUN_RADIUS_ENABLED=false
        ;;
    *)
        echo "ERROR: Invalid --reduced_based_run_radius enable value '${REDUCED_BASED_RUN_RADIUS_ENABLED_RAW}'."
        exit 1
        ;;
esac

if ${REDUCED_BASED_RUN_RADIUS_ENABLED}; then
    if ! [[ "${REDUCED_BASED_RUN_RADIUS_RAW}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
        echo "ERROR: --reduced_based_run_radius must be a positive number (meters). Got '${REDUCED_BASED_RUN_RADIUS_RAW}'"
        exit 1
    fi
    if ! awk "BEGIN {exit !(${REDUCED_BASED_RUN_RADIUS_RAW} > 0)}"; then
        echo "ERROR: --reduced_based_run_radius must be > 0. Got '${REDUCED_BASED_RUN_RADIUS_RAW}'"
        exit 1
    fi
    REDUCED_BASED_RUN_RADIUS="${REDUCED_BASED_RUN_RADIUS_RAW}"
else
    REDUCED_BASED_RUN_RADIUS="15"
fi

# =============================================================================
# SECTION 1: DERIVED PATHS AND NAMES (TREE STRUCTURE)
# =============================================================================

# Tree levels
TREE_ROOT="pdg${PDG}_E${ENERGY}_r${TEL_RADIUS}_s${SEED}"
TREE_ZENITH="zen${ZENITH}"
TREE_AZIMUTH="az${AZIMUTH}"
TREE_HEIGHT="h${INJ_HEIGHT}"
TREE_POSITION="x${TEL_X}_y${TEL_Y}_z${TEL_Z}"

# Full path to this run's directory
SIM_DIR="${OUTPUT_BASE_DIR}/${TREE_ROOT}/${TREE_ZENITH}/${TREE_AZIMUTH}/${TREE_HEIGHT}/${TREE_POSITION}"

# Base run path (x=0, y=0, z=0) — same tree but different leaf
BASE_POSITION="x0_y0_z0"
BASE_SIM_DIR="${OUTPUT_BASE_DIR}/${TREE_ROOT}/${TREE_ZENITH}/${TREE_AZIMUTH}/${TREE_HEIGHT}/${BASE_POSITION}"

# Flat tag used for the CHASM subdirectory name (container expects this)
GEOM_TAG="Tilt_pdg${PDG}_E${ENERGY}_zen${ZENITH}_az${AZIMUTH}_h${INJ_HEIGHT}_x${TEL_X}_y${TEL_Y}_z${TEL_Z}_r${TEL_RADIUS}_s${SEED}"
BASE_GEOM_TAG="Tilt_pdg${PDG}_E${ENERGY}_zen${ZENITH}_az${AZIMUTH}_h${INJ_HEIGHT}_x0_y0_z0_r${TEL_RADIUS}_s${SEED}"

# The base .dat file lives inside the base run's CORSIKA8 directory
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

# CSV path for combo completion status (used to force rerun when file_found=0)
CSV_STATUS_PATH="${OUTPUT_BASE_DIR}/csv_output/scan_care_pid${PDG}_E${ENERGY}_R${TEL_RADIUS}_y${TEL_Y}_s${SEED}.csv"

FORCE_RERUN_FROM_CSV=false
FORCE_RERUN_REASON=""
EFFECTIVE_RERUN_EVENT="${RERUN_EVENT}"

detect_force_rerun_from_csv() {
    if [ ! -f "${CSV_STATUS_PATH}" ]; then
        return
    fi

    local csv_status
    csv_status=$(python3 - "${CSV_STATUS_PATH}" "${SEED}" "${ZENITH}" "${AZIMUTH}" "${INJ_HEIGHT}" "${TEL_X}" "${TEL_Z}" <<'PYEOF'
import csv
import sys

path, seed, zen, az, height, tel_x, tel_z = sys.argv[1:]

def is_match(row):
    try:
        return (
            float(row["seed"]) == float(seed)
            and float(row["zen"]) == float(zen)
            and float(row["az"]) == float(az)
            and float(row["height"]) == float(height)
            and float(row["tel_x"]) == float(tel_x)
            and float(row["tel_z"]) == float(tel_z)
        )
    except Exception:
        return False

try:
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        required = {"seed", "zen", "az", "height", "tel_x", "tel_z", "file_found"}
        if not required.issubset(set(reader.fieldnames or [])):
            print("MISSING_HEADER")
            raise SystemExit

        for row in reader:
            if is_match(row):
                print((row.get("file_found", "") or "").strip())
                raise SystemExit

    print("NO_MATCH")
except Exception:
    print("CSV_READ_ERROR")
PYEOF
)

    case "${csv_status}" in
        0)
            FORCE_RERUN_FROM_CSV=true
            FORCE_RERUN_REASON="CSV combo has file_found=0 (${CSV_STATUS_PATH})"
            ;;
        1|NO_MATCH)
            ;;
        MISSING_HEADER)
            echo "WARNING: CSV status check skipped due to missing headers: ${CSV_STATUS_PATH}"
            ;;
        CSV_READ_ERROR)
            echo "WARNING: CSV status check failed to read ${CSV_STATUS_PATH}; continuing without force-rerun override."
            ;;
        *)
            echo "WARNING: CSV status check returned '${csv_status}' for ${CSV_STATUS_PATH}; continuing without force-rerun override."
            ;;
    esac
}

detect_force_rerun_from_csv
if ${FORCE_RERUN_FROM_CSV}; then
    EFFECTIVE_RERUN_EVENT=true
fi

# =============================================================================
# SECTION 1a: DETERMINE RUN MODE
# =============================================================================

IS_BASE_RUN=false
if [ "$(echo "${TEL_X} == 0 && ${TEL_Z} == 0" | bc -l)" -eq 1 ]; then
    IS_BASE_RUN=true
fi

echo ""
echo "============================================================"
echo "CORSIKA 8 + Trinity Chain (TREE LAYOUT)"
echo "SLURM Job ID       : ${SLURM_JOB_ID}"
echo "Tree path          : ${TREE_ROOT}/${TREE_ZENITH}/${TREE_AZIMUTH}/${TREE_HEIGHT}/${TREE_POSITION}"
echo "Geometry tag       : ${GEOM_TAG}"
echo "Base geometry tag  : ${BASE_GEOM_TAG}"
echo "SIM_DIR            : ${SIM_DIR}"
echo "BASE_SIM_DIR       : ${BASE_SIM_DIR}"
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
echo "Rerun event        : requested=${RERUN_EVENT}, effective=${EFFECTIVE_RERUN_EVENT}"
if ${FORCE_RERUN_FROM_CSV}; then
    echo "Rerun override     : ${FORCE_RERUN_REASON}"
fi
echo "Update timestamp   : ${UPDATE_TIME_STAMP}"
echo "Reduced base radius: enabled=${REDUCED_BASED_RUN_RADIUS_ENABLED}, value=${REDUCED_BASED_RUN_RADIUS}"
echo "Correction script  : ${CORRECTION_SCRIPT}"
if ${EFFECTIVE_RERUN_EVENT}; then
    echo "Run mode           : RERUN DOWNSTREAM (reuse cherenkov_hits_base.dat, no CORSIKA8)"
elif ${IS_BASE_RUN}; then
    echo "Run mode           : FULL RUN (base x=0, z=0)"
else
    echo "Run mode           : REUSE BASE (offset x=${TEL_X}, z=${TEL_Z})"
fi
echo "Is base run (x=0,z=0): ${IS_BASE_RUN}"
echo "Started at         : $(date)"
echo "============================================================"

# =============================================================================
# SECTION 1b: WRITE METADATA
# =============================================================================

# Ensure the SIM_DIR exists for metadata (will be created fully in Section 2)
mkdir -p "${SIM_DIR}"

cat > "${SIM_DIR}/metadata.yaml" <<EOF
pdg: ${PDG}
energy_eV: ${ENERGY}
zenith_deg: ${ZENITH}
azimuth_deg: ${AZIMUTH}
injection_height_m: ${INJ_HEIGHT}
observation_level_m: ${OBS_LEVEL}
telescope_x_m: ${TEL_X}
telescope_y_m: ${TEL_Y}
telescope_z_m: ${TEL_Z}
cherenkov_radius_m: ${CHERENKOV_RADIUS}
telescope_radius_m: ${TEL_RADIUS}
seed: ${SEED}
hadron_model: ${HADRON_MODEL}
rerun_event: ${EFFECTIVE_RERUN_EVENT}
rerun_event_requested: ${RERUN_EVENT}
update_time_stamp: ${UPDATE_TIME_STAMP}
is_base_run: ${IS_BASE_RUN}
tree_path: ${TREE_ROOT}/${TREE_ZENITH}/${TREE_AZIMUTH}/${TREE_HEIGHT}/${TREE_POSITION}
geom_tag: ${GEOM_TAG}
slurm_job_id: ${SLURM_JOB_ID}
created: $(date -Iseconds)
EOF

echo "Wrote metadata to: ${SIM_DIR}/metadata.yaml"

# =============================================================================
# SECTION 1c: CHECK IF BASE RUN ALREADY EXISTS (for x=0, z=0 only)
# =============================================================================

if ${IS_BASE_RUN} && ! ${EFFECTIVE_RERUN_EVENT}; then
    if ${UPDATE_TIME_STAMP}; then
        BASE_DAT_FOR_TOUCH=""
        for CANDIDATE in \
            "${BASE_DAT_CORSIKA8}" \
            "${BASE_DAT_CHASM}" \
            "${BASE_SIM_DIR}/CHASM/cherenkov_hits_base.dat" \
            "${BASE_SIM_DIR}/CORSIKA8/${BASE_GEOM_TAG}/cherenkov_hits_base.dat"; do
            if [ -f "${CANDIDATE}" ]; then
                BASE_DAT_FOR_TOUCH="${CANDIDATE}"
                break
            fi
        done

        if [ -n "${BASE_DAT_FOR_TOUCH}" ]; then
            if touch "${BASE_DAT_FOR_TOUCH}"; then
                echo "Updated timestamp: ${BASE_DAT_FOR_TOUCH}"
            else
                echo "WARNING: Failed to update timestamp: ${BASE_DAT_FOR_TOUCH}"
            fi
        else
            echo "Timestamp update requested, but no base cherenkov_hits_base.dat found yet."
        fi
    fi

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
# SECTION 1d: CHECK IF REUSE RUN ALREADY EXISTS
# =============================================================================

if ! ${IS_BASE_RUN} && ! ${EFFECTIVE_RERUN_EVENT}; then
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
# SECTION 1e: CLEAN UP PREVIOUS RUN OUTPUT (if any)
# =============================================================================

echo ""
echo "--- SECTION 1e: Cleaning up previous run output ---"
echo ""

if [ -d "${SIM_DIR}" ]; then
    echo "Previous SIM_DIR exists: ${SIM_DIR}"

    for OLD_DIR in \
        "${SIM_DIR}/CHASM" \
        "${SIM_DIR}/CARE" \
        "${SIM_DIR}/GROPT" \
        "${SIM_DIR}/CIO" \
        "${SIM_DIR}/corsika8_output"; do
        if [ -d "${OLD_DIR}" ]; then
            rm -rf "${OLD_DIR}"
            echo "  Removed: ${OLD_DIR}"
        fi
    done

    if ! ${EFFECTIVE_RERUN_EVENT} || ! ${IS_BASE_RUN}; then
        if [ -d "${SIM_DIR}/CORSIKA8" ]; then
            rm -rf "${SIM_DIR}/CORSIKA8"
            echo "  Removed: ${SIM_DIR}/CORSIKA8"
        fi
    fi

    for OLD_FILE in \
        "${SIM_DIR}/corsika8_tables.dat" \
        "${SIM_DIR}/cherenkov_hits_base.dat.preserved" \
        "${SIM_DIR}/CORSIKA8/cherenkov_hits.dat.orig" \
        "${SIM_DIR}/CHASM/cherenkov_hits.dat.orig"; do
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

    # Look for the base dat file in possible locations (tree structure)
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
        echo "  Base SIM_DIR: ${BASE_SIM_DIR}"
        echo ""
        echo "  The base run (TEL_X=0, TEL_Z=0) must complete first."
        exit 1
    fi

    echo "Found base dat file: ${BASE_DAT_FOUND}"
    echo "  Size: $(du -sh "${BASE_DAT_FOUND}" | cut -f1)"

    cp "${BASE_DAT_FOUND}" "${CHASM_TILT_DIR}/cherenkov_hits.dat"
    echo "Copied to: ${CHASM_TILT_DIR}/cherenkov_hits.dat"

else
    # -----------------------------------------------------------------
    # FULL RUN MODE: TEL_X == 0 AND TEL_Z == 0 (base run)
    # -----------------------------------------------------------------
    echo ""
    if ${EFFECTIVE_RERUN_EVENT}; then
        echo ""
        echo "--- SECTION 3: RERUN MODE (BASE) — Reusing cherenkov_hits_base.dat ---"
        echo ""

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
            echo "WARNING: effective rerun_event=true but cherenkov_hits_base.dat not found for base geometry."
            echo "  Checked: ${BASE_DAT_CORSIKA8}"
            echo "  Checked: ${BASE_DAT_CHASM}"
            echo "  Falling back to FULL RUN (run CORSIKA 8) for this base event."
            EFFECTIVE_RERUN_EVENT=false
        else
            echo "Found base dat file: ${BASE_DAT_FOUND}"
            echo "  Size: $(du -sh "${BASE_DAT_FOUND}" | cut -f1)"

            cp "${BASE_DAT_FOUND}" "${CHASM_TILT_DIR}/cherenkov_hits_base.dat"
            cp "${BASE_DAT_FOUND}" "${CHASM_TILT_DIR}/cherenkov_hits.dat"
            echo "Copied base dat to: ${CHASM_TILT_DIR}/cherenkov_hits_base.dat"
            echo "Copied working dat to: ${CHASM_TILT_DIR}/cherenkov_hits.dat"
        fi
    fi

    if ! ${EFFECTIVE_RERUN_EVENT}; then
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

        cp "${C8_DAT_FILE}" "${CHASM_TILT_DIR}/cherenkov_hits_base.dat"
        echo "Saved base dat:  ${CHASM_TILT_DIR}/cherenkov_hits_base.dat"

        mv "${C8_DAT_FILE}" "${CHASM_TILT_DIR}/cherenkov_hits.dat"
        echo "Moved dat to:    ${CHASM_TILT_DIR}/cherenkov_hits.dat"

        C8_OUTPUT_SUBDIR="${C8_WORK_DIR}/${C8_OUT_PREFIX}"
        if [ -d "${C8_OUTPUT_SUBDIR}" ]; then
            if [ -d "${SIM_DIR}/corsika8_output" ]; then
                rm -rf "${SIM_DIR}/corsika8_output"
            fi
            mv "${C8_OUTPUT_SUBDIR}" "${SIM_DIR}/corsika8_output"
            echo "CORSIKA 8 output dir moved to: ${SIM_DIR}/corsika8_output"
        fi

        if [ -f "${C8_WORK_DIR}/tables.dat" ]; then
            mv "${C8_WORK_DIR}/tables.dat" "${SIM_DIR}/corsika8_tables.dat"
            echo "Moved tables.dat to: ${SIM_DIR}/corsika8_tables.dat"
        fi

        cd "${SIM_DIR}"
        rm -rf "${C8_WORK_DIR}"
        echo "Removed work dir: ${C8_WORK_DIR}"
    fi
fi

echo ""
echo "CHASM tilt directory now contains:"
ls -lh "${CHASM_TILT_DIR}/"

# =============================================================================
# SECTION 4: PARTICLE LOCATION CORRECTION
# =============================================================================

echo ""
echo "--- SECTION 4: Particle location correction ---"
echo ""

CHERENKOV_DAT="${CHASM_TILT_DIR}/cherenkov_hits.dat"
CORRECTION_REPORT_JSON="${SIM_DIR}/correction_report_firstpass.json"
FIRSTPASS_REPORT_PREEXISTING=false
if ${IS_BASE_RUN} && [ -f "${CORRECTION_REPORT_JSON}" ]; then
    FIRSTPASS_REPORT_PREEXISTING=true
fi

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

CHERENKOV_BACKUP_TMP="${CHERENKOV_DAT}.tmp_backup"
cp "${CHERENKOV_DAT}" "${CHERENKOV_BACKUP_TMP}"
echo "Created temporary backup: ${CHERENKOV_BACKUP_TMP}"

python3 "${CORRECTION_SCRIPT}" \
    --cherenkov-eventio "${CHERENKOV_DAT}" \
    --telescope-radius "${TEL_RADIUS}" \
    --telescope-x "${TEL_X}" \
    --telescope-y "${TEL_Z}" \
    --output "${CHERENKOV_DAT}.corrected" \
    --recenter \
    --report-json "${CORRECTION_REPORT_JSON}"

CORR_EXIT=$?

if [ ${CORR_EXIT} -ne 0 ]; then
    echo "ERROR: Particle location correction failed with exit code ${CORR_EXIT}"
    echo "Restoring original file."
    mv "${CHERENKOV_BACKUP_TMP}" "${CHERENKOV_DAT}"
    exit ${CORR_EXIT}
fi

mv "${CHERENKOV_DAT}.corrected" "${CHERENKOV_DAT}"
rm -f "${CHERENKOV_BACKUP_TMP}"
echo ""
echo "Corrected file written: ${CHERENKOV_DAT}  ($(du -sh "${CHERENKOV_DAT}" | cut -f1))"
echo "Removed temporary backup: ${CHERENKOV_BACKUP_TMP}"

# Base-event metadata: store first-pass recenter location in telescope coordinates.
if ${IS_BASE_RUN}; then
    FIRST_PASS_RECENTER_X="${TEL_X}"
    FIRST_PASS_RECENTER_Y="${TEL_Z}"
    if [ -f "${CORRECTION_REPORT_JSON}" ]; then
        read -r FIRST_PASS_RECENTER_X FIRST_PASS_RECENTER_Y < <(
            python3 - "${CORRECTION_REPORT_JSON}" <<'PYEOF'
import json
import sys

with open(sys.argv[1], "r", encoding="utf-8") as f:
    d = json.load(f)

print(f"{d.get('center_x_m', 0.0)} {d.get('center_y_m', 0.0)}")
PYEOF
        )
    else
        echo "WARNING: correction report JSON not found, falling back to TEL_X/TEL_Z for metadata fields."
    fi
    {
        echo "recentered_telescope_x_m: ${FIRST_PASS_RECENTER_X}"
        echo "recentered_telescope_y_m: ${FIRST_PASS_RECENTER_Y}"
    } >> "${SIM_DIR}/metadata.yaml"
    echo "Updated base metadata with first-pass recenter position: x=${FIRST_PASS_RECENTER_X}, y=${FIRST_PASS_RECENTER_Y}"
fi

if ${REDUCED_BASED_RUN_RADIUS_ENABLED}; then
    if [ ! -f "${CORRECTION_REPORT_JSON}" ]; then
        echo "WARNING: --reduced_based_run_radius requested but correction_report_firstpass.json is missing at ${CORRECTION_REPORT_JSON}. Continuing."
    fi

    if ! ${IS_BASE_RUN}; then
        echo ""
        echo "Reduced base radius requested, but this is not a base run (x=0,z=0). Skipping reduction."
    else
        BASE_DAT_TARGET="${CHASM_TILT_DIR}/cherenkov_hits_base.dat"
        BASE_DAT_BACKUP_TMP="${BASE_DAT_TARGET}.tmp_backup"

        SHOULD_SKIP_REDUCTION=false
        if ${FIRSTPASS_REPORT_PREEXISTING} && [ -f "${BASE_DAT_TARGET}" ] && { ${EFFECTIVE_RERUN_EVENT} || ${UPDATE_TIME_STAMP}; }; then
            SHOULD_SKIP_REDUCTION=true
        fi

        if ${SHOULD_SKIP_REDUCTION}; then
            echo ""
            echo "Reduced base radius requested, but existing correction_report_firstpass.json and base dat were found."
            echo "Assuming base event is already reduced; skipping base-dat reduction pass."
            if ${UPDATE_TIME_STAMP}; then
                if touch "${BASE_DAT_TARGET}"; then
                    echo "Updated timestamp for reused base dat: ${BASE_DAT_TARGET}"
                else
                    echo "WARNING: Failed to update timestamp for reused base dat: ${BASE_DAT_TARGET}"
                fi
                if touch "${CORRECTION_REPORT_JSON}"; then
                    echo "Updated timestamp for reused correction report: ${CORRECTION_REPORT_JSON}"
                else
                    echo "WARNING: Failed to update timestamp for reused correction report: ${CORRECTION_REPORT_JSON}"
                fi
            fi
        else
        if ! awk "BEGIN {exit !(${REDUCED_BASED_RUN_RADIUS} <= ${CHERENKOV_RADIUS})}"; then
            echo "ERROR: --reduced_based_run_radius (${REDUCED_BASED_RUN_RADIUS} m) must be <= original base run radius (${CHERENKOV_RADIUS} m)."
            exit 1
        fi

        if [ ! -f "${BASE_DAT_TARGET}" ]; then
            echo "ERROR: reduced base radius requested, but base dat not found: ${BASE_DAT_TARGET}"
            exit 1
        fi

        echo ""
        echo "--- SECTION 4a: Reducing base cherenkov_hits_base.dat radius ---"
        echo "Base dat input:      ${BASE_DAT_TARGET}  ($(du -sh "${BASE_DAT_TARGET}" | cut -f1))"
        echo "Original radius:     ${CHERENKOV_RADIUS} m"
        echo "Reduced base radius: ${REDUCED_BASED_RUN_RADIUS} m"

        cp "${BASE_DAT_TARGET}" "${BASE_DAT_BACKUP_TMP}"

        python3 "${CORRECTION_SCRIPT}" \
            --cherenkov-eventio "${BASE_DAT_TARGET}" \
            --telescope-radius "${REDUCED_BASED_RUN_RADIUS}" \
            --telescope-x 0 \
            --telescope-y 0 \
            --output "${BASE_DAT_TARGET}.reduced"

        REDUCE_EXIT=$?
        if [ ${REDUCE_EXIT} -ne 0 ]; then
            echo "ERROR: Base dat reduction failed with exit code ${REDUCE_EXIT}"
            echo "Restoring original base dat."
            mv "${BASE_DAT_BACKUP_TMP}" "${BASE_DAT_TARGET}"
            exit ${REDUCE_EXIT}
        fi

        mv "${BASE_DAT_TARGET}.reduced" "${BASE_DAT_TARGET}"
        rm -f "${BASE_DAT_BACKUP_TMP}"

        echo "Reduced base dat written: ${BASE_DAT_TARGET}  ($(du -sh "${BASE_DAT_TARGET}" | cut -f1))"
        fi
    fi
fi

# =============================================================================
# SECTION 5: PERMISSIONS BEFORE CONTAINER
# =============================================================================

echo ""
echo "--- SECTION 5: Setting permissions ---"

BASE_DAT_IN_CHASM="${CHASM_TILT_DIR}/cherenkov_hits_base.dat"
BASE_DAT_TEMP="${SIM_DIR}/cherenkov_hits_base.dat.tmp"

if [ -f "${BASE_DAT_IN_CHASM}" ]; then
    mv "${BASE_DAT_IN_CHASM}" "${BASE_DAT_TEMP}"
    echo "Temporarily moved cherenkov_hits_base.dat out of CHASM tilt dir"
    echo "  From: ${BASE_DAT_IN_CHASM}"
    echo "  To:   ${BASE_DAT_TEMP}"
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

if [ -f "${BASE_DAT_TEMP}" ]; then
    mv "${BASE_DAT_TEMP}" "${BASE_DAT_IN_CHASM}"
    echo "Restored cherenkov_hits_base.dat to: ${BASE_DAT_IN_CHASM}"
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
    if [ -d "${SIM_DIR}/CORSIKA8" ]; then
        rm -rf "${SIM_DIR}/CORSIKA8"
        echo "Removed existing ${SIM_DIR}/CORSIKA8 before rename"
    fi
    mv "${SIM_DIR}/CHASM" "${SIM_DIR}/CORSIKA8"
    echo "Renamed: ${SIM_DIR}/CHASM → ${SIM_DIR}/CORSIKA8"
else
    echo "WARNING: ${SIM_DIR}/CHASM not found, skipping rename."
fi

echo ""
echo "--- Removing backup-copy artifacts (.orig/.preserved) ---"
find "${SIM_DIR}" -type f \( \
    -name "cherenkov_hits.dat.orig" -o \
    -name "cherenkov_hits_base.dat.preserved" -o \
    -name "cherenkov_hits.dat.tmp_backup" -o \
    -name "cherenkov_hits_base.dat.tmp_backup" \
\) -print -delete || true

# =============================================================================
# SECTION 8: FINAL STATUS
# =============================================================================

cd "${SIM_DIR}"

echo ""
echo "============================================================"
if [ ${CONTAINER_EXIT} -eq 0 ]; then
    echo "SUCCESS: Full chain completed."
    echo ""
    echo "Final tree location: ${TREE_ROOT}/${TREE_ZENITH}/${TREE_AZIMUTH}/${TREE_HEIGHT}/${TREE_POSITION}"
    echo ""
    echo "Final output structure:"
    find "${SIM_DIR}" -maxdepth 4 -type f | sort
else
    echo "ERROR: Container chain failed with exit code ${CONTAINER_EXIT}."
fi
echo "Finished at: $(date)"
echo "============================================================"

exit ${CONTAINER_EXIT}
