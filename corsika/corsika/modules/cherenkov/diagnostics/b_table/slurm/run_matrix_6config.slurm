#!/bin/bash
# =============================================================================
# run_matrix_6config.slurm
#
# The six-configuration Cherenkov matrix (2g-h2n, upgoing muon), one SLURM job
# per configuration via a job array.
#
#   0 straight_noabs        plain geometry, no absorption
#   1 straight_abs          plain geometry + tabulated absorption (chord tau)
#   2 curve_noabs           curved ray, geometry table, no absorption
#   3 curve_abs             curved ray + geometry table + ABSORPTION TABLE
#   4 curve_noabs_notable   curved ray via the legacy tracer, no absorption
#   5 curve_abs_notable     curved ray via the legacy tracer + path-integral tau
#
# Physics settings default to the documented 2g-h2n baseline (zenith 92,
# telescope radius 15 m, MODTRAN_new absorption file); both are overridable so
# the same harness can reproduce any configuration.  Only the three flags that
# select the route differ between configurations, so the photon stream is
# IDENTICAL across all six (seed 1000, absorption is deterministic and draws no
# random numbers) - which is what makes a per-photon A/B comparison meaningful.
#
# The table cache is shared across zenith/radius: the geometry and absorption
# tables depend on the observation level and the atmosphere profile, NOT on where
# the shower points or how big the dish is.  So every configuration and every
# energy hits the same cached build.
#
# Usage:  sbatch run_matrix_6config.slurm [ENERGY] [ZENITH] [TEL_RADIUS_M]
# =============================================================================
#SBATCH --account=owner-guest
#SBATCH --partition=kingspeak-guest
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --job-name=cher_matrix
#SBATCH --array=0-5
#SBATCH --output=/scratch/general/vast/u1520754/log/cher_matrix_%A_%a.out
#SBATCH --error=/scratch/general/vast/u1520754/log/cher_matrix_%A_%a.err

set -uo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
APP=/uufs/chpc.utah.edu/common/home/u1520754/corsika-build/applications/muon_cherenkov_3d
ABS=/uufs/chpc.utah.edu/common/home/u1520754/corsika/modules/data/CHERENKOV/atmosphere/abstable_MODTRAN_new.dat
ENERGY="${1:-1e5}"
ZENITH="${2:-92}"
TEL_R="${3:-15}"
BASE=/scratch/general/vast/u1520754/cher_wave5/zen${ZENITH}_r${TEL_R}_E${ENERGY}
CACHE=/scratch/general/vast/u1520754/cher_table_cache/obs3000

CONFIGS=(straight_noabs straight_abs curve_noabs curve_abs curve_noabs_notable curve_abs_notable)
NAME="${CONFIGS[$SLURM_ARRAY_TASK_ID]}"
OUT="${BASE}/${NAME}"

mkdir -p "${BASE}" "${CACHE}" /scratch/general/vast/u1520754/log
rm -rf "${OUT}"
mkdir -p "${OUT}"

# per-configuration route flags (everything else is identical)
EXTRA=()
case "${NAME}" in
  straight_noabs)      ;;                                   # no curvature, no absorption
  straight_abs)        EXTRA+=(--atmosphere-absorption-file "${ABS}") ;;
  curve_noabs)         EXTRA+=(--enable-curvature) ;;
  curve_abs)           EXTRA+=(--enable-curvature --atmosphere-absorption-file "${ABS}") ;;
  curve_noabs_notable) EXTRA+=(--enable-curvature --disable-cherenkov-curved-table) ;;
  curve_abs_notable)   EXTRA+=(--enable-curvature --atmosphere-absorption-file "${ABS}" \
                               --disable-cherenkov-curved-table) ;;
esac

echo "======================================================================"
echo " config   : ${NAME}   (array task ${SLURM_ARRAY_TASK_ID})"
echo " energy   : ${ENERGY}   zenith ${ZENITH}   radius ${TEL_R} m"
echo " host     : $(hostname)   job ${SLURM_JOB_ID}"
echo " output   : ${OUT}"
echo " flags    : ${EXTRA[*]:-none}"
echo " started  : $(date)"
echo "======================================================================"

START=$(date +%s)
"${APP}" \
  --pdg 15 -E "${ENERGY}" --injection-height 1500 \
  --emcut 0.015 --mucut 0.015 --hadcut 0.3 \
  --observation-level 3000 --zenith "${ZENITH}" --azimuth 270 \
  --cherenkov-projection telescope \
  --telescope-azimuth 270 --telescope-zenith 90 \
  --telescope-x 0 --telescope-y 0 --telescope-z 0 --telescope-radius "${TEL_R}" \
  --seed 1000 --hadronModel SIBYLL-2.3d \
  --cherenkov-output-format eventio \
  --cherenkov-table-cache "${CACHE}" \
  -f "${OUT}/out" -v info \
  "${EXTRA[@]}"
RC=$?
WALL=$(( $(date +%s) - START ))

echo ""
echo "--- key diagnostics (${NAME}) ---"
echo "WALL_SECONDS=${WALL}  APP_RC=${RC}  CONFIG=${NAME}  ENERGY=${ENERGY}"
exit ${RC}
