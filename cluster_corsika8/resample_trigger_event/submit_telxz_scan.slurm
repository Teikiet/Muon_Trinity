#!/bin/bash
#SBATCH --account=owner-guest
#SBATCH --partition=kingspeak-guest
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=3G
#SBATCH --job-name=telxz_submit
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1520754/logs/telxz_submit_%j.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/u1520754/logs/telxz_submit_%j.err

set -euo pipefail

mkdir -p /scratch/general/vast/u1520754/log

bash $HOME/Muon_Trinity/cluster_corsika8/resample_trigger_event/submit_triggered_telxz_scan.sh \
  --tel-x-values "0" \
  --tel-z-values "0" \
  --pe-threshold 20 \
  --seeds "1" \
  --rerun-event
