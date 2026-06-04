#!/bin/bash
#SBATCH --job-name=MCEq
#SBATCH --time=70:00:00
#SBATCH --ntasks=1
#SBATCH --mem=6GB

# Load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate jupyter_env
python /uufs/chpc.utah.edu/common/home/u1520754/Muon_Trinity/analysis/make_MCEq_grid.py


