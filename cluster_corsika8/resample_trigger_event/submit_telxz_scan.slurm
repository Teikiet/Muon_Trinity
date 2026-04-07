#!/bin/bash
#SBATCH --account=owner-guest
#SBATCH --partition=kingspeak-guest
#SBATCH --time=0:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --job-name=telxz_submit
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1520754/logs/telxz_submit_%j.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/u1520754/logs/telxz_submit_%j.err

set -euo pipefail

mkdir -p /scratch/general/vast/u1520754/log

bash $HOME/Muon_Trinity/cluster_corsika8/resample_trigger_event/submit_triggered_telxz_scan.sh \
  --tel-x-values "-5 -2 -1 -0.3 0 0.3 1 2 5" \
  --tel-z-values "-5 -2 -1 -0.3 0 0.3 1 2 5" \
  --pe-threshold 20 \
  --seeds "1"
