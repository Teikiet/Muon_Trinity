#!/bin/bash
#SBATCH --job-name=MCEq
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16GB
#SBATCH --partition=kingspeak-guest
#SBATCH --account=owner-guest
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1520754/logs/MCEq.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/u1520754/logs/MCEq.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate jupyter_env

python -u /uufs/chpc.utah.edu/common/home/u1520754/Muon_Trinity/analysis/make_MCEq_grid.py \
    > /uufs/chpc.utah.edu/common/home/u1520754/logs/MCEq.log 2>&1
