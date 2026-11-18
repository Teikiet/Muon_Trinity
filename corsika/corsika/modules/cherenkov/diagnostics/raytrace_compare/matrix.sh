#!/bin/bash
# matrix.sh -- six ray-tracing configurations x 5 seeds, the user's verified muon point.
# Diagnostics only; writes ONLY under /scratch/general/vast/u1520754/chercmp/
#   A straight line, no absorption (baseline)
#   B straight line + absorption
#   C curved legacy tracer, no absorption   (forces the tracer, bypasses the table)
#   D curved legacy tracer + absorption
#   E curved + CurvedRayTable, no absorption
#   F straight line, no absorption, LEGACY segment-loop timing
set -u
BIN=/uufs/chpc.utah.edu/common/home/u1520754/corsika-build/applications/muon_cherenkov_3d
ABS=/uufs/chpc.utah.edu/common/home/u1520754/corsika/modules/data/CHERENKOV/atmosphere/abstable_MODTRAN_new.dat
OUT=/scratch/general/vast/u1520754/chercmp
mkdir -p "$OUT"
cd "$OUT"

BASE=(--pdg 13 -E 1e4 --injection-height 43000 --emcut 0.02 --mucut 0.43 --hadcut 0.5
      --observation-level 2944 --zenith 90 --azimuth 270
      --cherenkov-projection telescope --telescope-azimuth 270 --telescope-zenith 90
      --telescope-x 0 --telescope-y 0 --telescope-z 0 --telescope-radius 15.0
      --telescope-position-correction 1
      --hadronModel SIBYLL-2.3d --cherenkov-output-format eventio --cherenkov-profile -v info)

run() {  # run <name> <seed> [extra flags...]
  local name=$1 seed=$2; shift 2
  "$BIN" "${BASE[@]}" -s "$seed" "$@" -f "$OUT/m_${name}_s${seed}" \
      > "$OUT/m_${name}_s${seed}.log" 2>&1
  printf '  %-2s seed=%-4s exit=%s hits=%s\n' "$name" "$seed" "$?" \
      "$(grep -o 'Total Cherenkov hits recorded: [0-9]*' "$OUT/m_${name}_s${seed}.log" | grep -o '[0-9]*$')"
}

for seed in 600 601 602 603 604; do
  echo "== seed $seed"
  run A "$seed"
  run B "$seed" --atmosphere-absorption-file "$ABS"
  run C "$seed" --enable-curvature --disable-cherenkov-curved-table
  run D "$seed" --enable-curvature --atmosphere-absorption-file "$ABS"
  run E "$seed" --enable-curvature
  run F "$seed" --disable-cherenkov-fast-timing
done
echo "MATRIX COMPLETE"
