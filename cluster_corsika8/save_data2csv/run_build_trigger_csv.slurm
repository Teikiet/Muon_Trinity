#!/bin/bash
#SBATCH --account=owner-guest
#SBATCH --partition=kingspeak-guest
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --job-name=build_trigger_csv
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1520754/csv_logs/build_trigger_csv_%j.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/u1520754/csv_logs/build_trigger_csv_%j.err

set -euo pipefail

BASE_PATH="/scratch/general/vast/u1520754/muon_sim_chain_tree"
SCRIPT="$HOME/Muon_Trinity/cluster_corsika8/save_data2csv/build_trigger_resample_csv.py"
PDG=13
RADIUS=5
TEL_Y=0
SEEDS=(1 2 3)
PE_THRESHOLD=20
TEL_X_VALUES="-5 -2 -1 -0.3 0 0.3 1 2 5"
TEL_Z_VALUES="-5 -2 -1 -0.3 0 0.3 1 2 5"

mkdir -p "$HOME/csv_logs"

# Generate energy strings
source /uufs/chpc.utah.edu/common/home/u1520754/miniconda3/etc/profile.d/conda.sh
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
E = E[(E >= 1e3) & (E <= 1e4)]

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

for ENERGY_STR in "${ENERGY_STRS[@]}"; do
    for SEED in "${SEEDS[@]}"; do
        echo "--- E=${ENERGY_STR} seed=${SEED} ---"

        python3 "$SCRIPT" \
            --base-path "$BASE_PATH" \
            --pid "$PDG" \
            --energy-str "$ENERGY_STR" \
            --radius "$RADIUS" \
            --tel-y "$TEL_Y" \
            --seed "$SEED" \
            --pe-threshold "$PE_THRESHOLD" \
            --tel-x-values $TEL_X_VALUES \
            --tel-z-values $TEL_Z_VALUES \
        || echo "  WARN: Failed for E=${ENERGY_STR} seed=${SEED}"

    done
done

echo "Done."
