#!/bin/bash
# =============================================================================
# run_clearance_sweep.slurm
#
# Sweeps the grazing rejection band (--cherenkov-grazing-clearance) and records the
# yields, so the IMAGE difference against the straight line can be computed as a
# function of the band:
#
#     net_noabs(c) = (curve_noabs(c) - straight_noabs) / straight_noabs
#     net_abs(c)   = (curve_abs(c)   - straight_abs)   / straight_abs
#
# The band is a LOOKUP guard, so it does not enter the table or its cache key --
# one cached build serves every value, which is why the sweep needs no rebuilds.
#
# Usage: sbatch run_clearance_sweep.slurm [ENERGY] [ZENITH] [TEL_RADIUS_M]
# =============================================================================
#SBATCH --account=owner-guest
#SBATCH --partition=kingspeak-guest
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --job-name=cher_sweep
#SBATCH --array=0-7
#SBATCH --output=/scratch/general/vast/u1520754/log/cher_sweep_%A_%a.out
#SBATCH --error=/scratch/general/vast/u1520754/log/cher_sweep_%A_%a.err

set -uo pipefail

APP=/uufs/chpc.utah.edu/common/home/u1520754/corsika-build/applications/muon_cherenkov_3d
ABS=/uufs/chpc.utah.edu/common/home/u1520754/corsika/modules/data/CHERENKOV/atmosphere/abstable_MODTRAN_new.dat
ENERGY="${1:-1e5}"
ZENITH="${2:-92}"
TEL_R="${3:-15}"
BASE=/scratch/general/vast/u1520754/cher_wave5/zen${ZENITH}_r${TEL_R}_E${ENERGY}
CACHE=/scratch/general/vast/u1520754/cher_table_cache/obs3000_v4

CLEAR=(10 100 200 400)
TASK=${SLURM_ARRAY_TASK_ID}
CI=$(( TASK / 2 ))
if [ $(( TASK % 2 )) -eq 0 ]; then NAME=curve_noabs; else NAME=curve_abs; fi
C=${CLEAR[$CI]}
OUT="${BASE}/sweep/${NAME}_c${C}"

mkdir -p "${BASE}/sweep" "${CACHE}" /scratch/general/vast/u1520754/log
rm -rf "${OUT}"
mkdir -p "${OUT}"

EXTRA=()
[ "${NAME}" = "curve_abs" ] && EXTRA+=(--atmosphere-absorption-file "${ABS}")

echo "=== ${NAME}  clearance=${C} m  energy=${ENERGY} on $(hostname) ==="
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
  --enable-curvature --cherenkov-grazing-clearance "${C}" \
  -f "${OUT}/out" -v info "${EXTRA[@]}"
RC=$?
echo "SWEEP name=${NAME} clearance=${C} WALL_SECONDS=$(( $(date +%s) - START )) APP_RC=$RC"
exit ${RC}
