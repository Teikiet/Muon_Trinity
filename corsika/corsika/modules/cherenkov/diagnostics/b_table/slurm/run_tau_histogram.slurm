#!/bin/bash
# =============================================================================
# run_tau_histogram.slurm
#
# Settles the W2 tau question.  Runs the two absorption routes back to back with
# the in-process tau/survival histogram enabled, so the printed line gives:
#   N_tau, tau_mean, tau_rms, tau_max  -- what reached the weight
#   surv_mean, implied absorbed %      -- what the weight actually did
#   surv buckets                       -- whether the distribution is skewed
#
# If the two routes' tau_means differ by a factor of ~2 -> the tau really is the
# difference and it lives in the distribution.  If the tau_means agree but the
# absorbed fractions differ -> the bookkeeping, not the physics.
#
# Usage:  sbatch run_tau_histogram.slurm [ENERGY] [ZENITH] [TEL_RADIUS_M]
# =============================================================================
#SBATCH --account=owner-guest
#SBATCH --partition=kingspeak-guest
#SBATCH --time=00:40:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --job-name=cher_tauhist
#SBATCH --array=0-1
#SBATCH --output=/scratch/general/vast/u1520754/log/cher_tauhist_%A_%a.out
#SBATCH --error=/scratch/general/vast/u1520754/log/cher_tauhist_%A_%a.err

set -uo pipefail

APP=/uufs/chpc.utah.edu/common/home/u1520754/corsika-build/applications/muon_cherenkov_3d
ABS=/uufs/chpc.utah.edu/common/home/u1520754/corsika/modules/data/CHERENKOV/atmosphere/abstable_MODTRAN_new.dat
ENERGY="${1:-1e5}"
ZENITH="${2:-92}"
TEL_R="${3:-15}"
BASE=/scratch/general/vast/u1520754/cher_wave5/zen${ZENITH}_r${TEL_R}_E${ENERGY}
CACHE=/scratch/general/vast/u1520754/cher_table_cache/obs3000

CONFIGS=(curve_abs curve_abs_notable)
NAME="${CONFIGS[$SLURM_ARRAY_TASK_ID]}"
OUT="${BASE}/${NAME}_tauhist"
mkdir -p "${BASE}" "${CACHE}" /scratch/general/vast/u1520754/log
rm -rf "${OUT}"
mkdir -p "${OUT}"

EXTRA=()
[ "${NAME}" = "curve_abs_notable" ] && EXTRA+=(--disable-cherenkov-curved-table)

echo "=== ${NAME} (energy ${ENERGY}, zenith ${ZENITH}, r ${TEL_R}) on $(hostname) ==="
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
  --enable-curvature --atmosphere-absorption-file "${ABS}" \
  -f "${OUT}/out" -v info "${EXTRA[@]}"
RC=$?
echo "WALL_SECONDS=$(( $(date +%s) - START ))  APP_RC=$RC  CONFIG=${NAME}"
exit ${RC}
