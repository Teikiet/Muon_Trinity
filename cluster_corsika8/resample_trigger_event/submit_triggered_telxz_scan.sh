#!/bin/bash
# =============================================================================
# submit_triggered_telxz_scan.sh
#
# Reads merged CSVs, filters for max_pe >= 20, and submits downstream-only
# jobs for a grid of (tel_x, tel_z) offsets.
#
# Usage:
#   bash submit_triggered_telxz_scan.sh \
#     --tel-x-values "-100 -50 0 50 100" \
#     --tel-z-values "-100 -50 0 50 100" \
#     [--pe-threshold 20] [--dry-run] [--rerun-event]
# =============================================================================

set -euo pipefail

# =============================================================================
# DEFAULTS
# =============================================================================
PE_THRESHOLD=20
DRY_RUN=false
RERUN_EVENT=false
TEL_X_VALUES=""
TEL_Z_VALUES=""

PDG=13
TEL_Y=0
TEL_RADIUS=5
OBS_LEVEL=2944
HADRON_MODEL="SIBYLL-2.3d"
SEEDS=(1)

BASE_PATH="/scratch/general/vast/u1520754/muon_sim_chain_tree"
SLURM_SCRIPT="$HOME/Muon_Trinity/cluster_corsika8/base_muon_event_generation/run_corsika8_trinity_chain_tree.slurm"

MAX_JOBS=500
SLEEP_SEC=30
MAX_SUBMIT_RETRIES=20

# =============================================================================
# ARGUMENT PARSING
# =============================================================================
while [[ $# -gt 0 ]]; do
    case "$1" in
        --tel-x-values)   TEL_X_VALUES="$2"; shift 2 ;;
        --tel-z-values)   TEL_Z_VALUES="$2"; shift 2 ;;
        --pe-threshold)   PE_THRESHOLD="$2"; shift 2 ;;
        --dry-run)        DRY_RUN=true; shift ;;
        --rerun-event)    RERUN_EVENT=true; shift ;;
        --seeds)          IFS=' ' read -ra SEEDS <<< "$2"; shift 2 ;;
        --max-jobs)       MAX_JOBS="$2"; shift 2 ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

if [ -z "$TEL_X_VALUES" ] || [ -z "$TEL_Z_VALUES" ]; then
    echo "ERROR: Must specify --tel-x-values and --tel-z-values"
    echo "Example: --tel-x-values '-100 -50 0 50 100' --tel-z-values '-100 -50 0 50 100'"
    exit 1
fi

IFS=' ' read -ra TEL_XS <<< "$TEL_X_VALUES"
IFS=' ' read -ra TEL_ZS <<< "$TEL_Z_VALUES"

echo "============================================================"
echo "Triggered tel_x/tel_z scan submission"
echo "============================================================"
echo "PE threshold:  ${PE_THRESHOLD}"
echo "tel_x values:  ${TEL_XS[*]}"
echo "tel_z values:  ${TEL_ZS[*]}"
echo "Seeds:         ${SEEDS[*]}"
echo "Dry run:       ${DRY_RUN}"
echo "Rerun event:   ${RERUN_EVENT}"
echo "============================================================"

# =============================================================================
# HELPER: compute cherenkov radius
# =============================================================================
compute_cherenkov_radius() {
    local ENERGY="$1"
    python3 -c "
import math
E = float('${ENERGY}')
r = max(15.0, 1.5e6 / E)
print(f'{r:.1f}')
"
}

# =============================================================================
# HELPER: submit with retry
# =============================================================================


submit_with_retry() {
    local retries=0
    local JID=""
    while [ $retries -lt $MAX_SUBMIT_RETRIES ]; do
        JID=$(sbatch --parsable "$@" 2>/tmp/sbatch_err_$$)
        local rc=$?
        local errmsg=$(cat /tmp/sbatch_err_$$ 2>/dev/null)
        rm -f /tmp/sbatch_err_$$

        if [ $rc -eq 0 ] && [ -n "$JID" ]; then
            echo "$JID"
            return 0
        fi

        if echo "$errmsg" | grep -qi "QOSMaxSubmitJobPerUserLimit\|Resource temporarily unavailable\|Socket timed out"; then
            retries=$((retries + 1))
            echo "  Retry $retries/$MAX_SUBMIT_RETRIES (${errmsg})..." >&2
            sleep $SLEEP_SEC
        else
            echo "  FATAL sbatch error: ${errmsg}" >&2
            return 1
        fi
    done
    echo "  Exhausted retries" >&2
    return 1
}

# =============================================================================
# HELPER: wait for job slots
# =============================================================================
wait_for_slots() {
    while true; do
        local n_running
        n_running=$(squeue -u "$USER" -h | wc -l)
        if [ "$n_running" -lt "$MAX_JOBS" ]; then
            return
        fi
        echo "  Queue full ($n_running >= $MAX_JOBS), waiting ${SLEEP_SEC}s..."
        sleep "$SLEEP_SEC"
    done
}

# =============================================================================
# STEP 1: Generate energy list from MCEq
# =============================================================================
echo ""
echo "--- Generating energy grid from MCEq ---"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate jupyter_env

ENERGY_FILE=$(mktemp)
python3 <<'PYEOF' > "$ENERGY_FILE"
import sys, os, math
os.environ["MCEQ_LOG_LEVEL"] = "50"
import logging
logging.disable(logging.CRITICAL)
import MCEq.config as config
config.kernel_config = 'MKL'
config.e_min = 1000
config.integrator = 'euler'

old_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')
from MCEq.core import MCEqRun
import crflux.models as pm
mceq = MCEqRun(
    interaction_model='SIBYLL23C',
    primary_model=(pm.GlobalSplineFitBeta, None),
    theta_deg=0.,
    density_model=("CORSIKA", ('USStd', None)),
)
sys.stdout = old_stdout

E = mceq.e_grid
E = E[(E > 1e3) & (E <= 1e4)]

def to_1e(val):
    exp = int(math.log10(val))
    coeff = val / 10**exp
    return f"{coeff:g}e{exp}"

for e in E:
    print(to_1e(e))
PYEOF

mapfile -t ENERGY_STRS < "$ENERGY_FILE"
rm -f "$ENERGY_FILE"
echo "Found ${#ENERGY_STRS[@]} energies"

# =============================================================================
# STEP 2: Scan all CSVs and extract triggered events
# =============================================================================
echo ""
echo "--- Extracting triggered events (max_pe >= ${PE_THRESHOLD}) ---"

TRIGGERED_FILE=$(mktemp)

python3 - "$BASE_PATH" "$PDG" "$TEL_RADIUS" "$TEL_Y" "$PE_THRESHOLD" \
    "${SEEDS[*]}" "${ENERGY_STRS[*]}" "$TRIGGERED_FILE" <<'PYEOF'
import sys, os, csv

base_path    = sys.argv[1]
pdg          = sys.argv[2]
tel_radius   = sys.argv[3]
tel_y        = sys.argv[4]
pe_threshold = float(sys.argv[5])
seeds        = sys.argv[6].split()
energies     = sys.argv[7].split()
output_file  = sys.argv[8]

triggered = []
total_checked = 0
total_triggered = 0

for energy_str in energies:
    for seed in seeds:
        csv_path = os.path.join(
            base_path,
            f"Muon_pid{pdg}_E{energy_str}_R{tel_radius}",
            "csv_output",
            f"scan_care_pid{pdg}_E{energy_str}_R{tel_radius}_y{tel_y}_s{seed}.csv"
        )
        if not os.path.exists(csv_path):
            continue

        with open(csv_path, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                total_checked += 1
                if row.get("file_found", "0") != "1":
                    continue

                max_pe = float(row.get("max_pe", 0))
                if max_pe < pe_threshold:
                    continue

                # Only select base events (tel_x=0, tel_z=0)
                if float(row.get("tel_x", 0)) != 0 or float(row.get("tel_z", 0)) != 0:
                    continue

                total_triggered += 1
                triggered.append({
                    "energy_string": row["energy_string"],
                    "seed":          row["seed"],
                    "zen":           row["zen"],
                    "az":            row["az"],
                    "height":        row["height"],
                    "max_pe":        row["max_pe"],
                })

print(f"Checked {total_checked} rows, found {total_triggered} triggered base events",
      file=sys.stderr)

# Write triggered list
with open(output_file, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["energy_string","seed","zen","az","height","max_pe"])
    writer.writeheader()
    for t in triggered:
        writer.writerow(t)

print(f"Wrote {len(triggered)} triggered events to {output_file}", file=sys.stderr)
PYEOF

N_TRIGGERED=$(tail -n +2 "$TRIGGERED_FILE" | wc -l)
echo "Triggered events: ${N_TRIGGERED}"

if [ "$N_TRIGGERED" -eq 0 ]; then
    echo "No triggered events found. Nothing to submit."
    rm -f "$TRIGGERED_FILE"
    exit 0
fi

# =============================================================================
# STEP 3: Submit downstream-only jobs for each triggered event × (tel_x, tel_z)
# =============================================================================
echo ""
echo "--- Submitting downstream-only jobs ---"

N_TEL_POSITIONS=$(( ${#TEL_XS[@]} * ${#TEL_ZS[@]} ))
N_NEW_POSITIONS=$((N_TEL_POSITIONS))
TOTAL_JOBS=$((N_TRIGGERED * N_NEW_POSITIONS))

echo "Telescope positions per event: ${N_NEW_POSITIONS}"
echo "Total jobs to submit: ${TOTAL_JOBS}"

SUBMITTED=0
SKIPPED=0
FAILED=0

# Shorter time limit for downstream-only jobs
DOWNSTREAM_TIME="1:00:00"

while IFS=, read -r energy_string seed zen az height max_pe; do
    # Skip header
    [ "$energy_string" = "energy_string" ] && continue

    # Fix: strip decimal from height (19090.0 -> 19090)
    height=$(printf '%.0f' "$height")

    # Energy group directory
    ENERGY_GROUP="Muon_pid${PDG}_E${energy_string}_R${TEL_RADIUS}"
    TREE_ROOT="pdg${PDG}_E${energy_string}_r${TEL_RADIUS}_s${seed}"

    # Convert energy string to eV for the slurm script
    ENERGY_GEV=$(python3 -c "print(float('${energy_string}'))")
    CHERENKOV_RADIUS=$(compute_cherenkov_radius "$ENERGY_GEV")

    for tel_x in "${TEL_XS[@]}"; do
        for tel_z in "${TEL_ZS[@]}"; do

            # Skip base position unless rerun mode explicitly requests it.
            if ! $RERUN_EVENT && [ "$tel_x" = "0" ] && [ "$tel_z" = "0" ]; then
                SKIPPED=$((SKIPPED + 1))
                continue
            fi

            # Check if output already exists
            TREE_POS="x${tel_x}_y${TEL_Y}_z${tel_z}"
            CHECK_DIR="${BASE_PATH}/${ENERGY_GROUP}/${TREE_ROOT}/zen${zen}/az${az}/h${height}/${TREE_POS}"

            if ! $RERUN_EVENT && [ -d "${CHECK_DIR}/CARE" ] && [ "$(find "${CHECK_DIR}/CARE" -name '*.root' 2>/dev/null | head -1)" ]; then
                SKIPPED=$((SKIPPED + 1))
                continue
            fi

            # Build log path
            LOG_DIR="${BASE_PATH}/Muon_pid${PDG}_E${energy_string}_R${TEL_RADIUS}/logs"
            mkdir -p "$LOG_DIR"
            LOG_FILE="${LOG_DIR}/downstream_pid${PDG}_E${energy_string}_zen${zen}_az${az}_h${height}_x${tel_x}_z${tel_z}_s${seed}.log"

            if $DRY_RUN; then
                echo "[DRY-RUN] Would submit: E=${energy_string} zen=${zen} az=${az} h=${height} x=${tel_x} z=${tel_z} s=${seed}"
                SUBMITTED=$((SUBMITTED + 1))
                continue
            fi
            wait_for_slots

            JID=$(submit_with_retry \
                --time="$DOWNSTREAM_TIME" \
                --mem=8G \
                "$SLURM_SCRIPT" \
                "${BASE_PATH}/${ENERGY_GROUP}" \
                "$LOG_FILE" \
                "$PDG" \
                "$energy_string" \
                "$zen" \
                "$az" \
                "$height" \
                "$OBS_LEVEL" \
                "$tel_x" \
                "$TEL_Y" \
                "$tel_z" \
                "$CHERENKOV_RADIUS" \
                "$TEL_RADIUS" \
                "$seed" \
                "$HADRON_MODEL" \
                "$RERUN_EVENT"
            ) || { FAILED=$((FAILED + 1)); continue; }

            SUBMITTED=$((SUBMITTED + 1))

            if (( SUBMITTED % 100 == 0 )); then
                echo "  Progress: submitted=${SUBMITTED} skipped=${SKIPPED} failed=${FAILED}"
            fi
        done
    done
done < "$TRIGGERED_FILE"

rm -f "$TRIGGERED_FILE"

echo ""
echo "============================================================"
echo "SUMMARY"
echo "============================================================"
echo "Triggered base events: ${N_TRIGGERED}"
echo "Tel positions per event: ${N_NEW_POSITIONS}"
echo "Submitted: ${SUBMITTED}"
echo "Skipped (exists or base unless --rerun-event): ${SKIPPED}"
echo "Failed: ${FAILED}"
echo "============================================================"
