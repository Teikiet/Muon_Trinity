#!/bin/bash
#SBATCH --job-name=care_cnn
#SBATCH --account=owner-guest
#SBATCH --partition=kingspeak-guest
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1520754/Muon_Trinity/logs/care_cnn_%j.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/u1520754/Muon_Trinity/logs/care_cnn_%j.err

# ── Info ──────────────────────────────────────────────────────────────────────
echo "Job started on $(hostname) at $(date)"
echo "SLURM_JOB_ID = $SLURM_JOB_ID"

mkdir -p /uufs/chpc.utah.edu/common/home/u1520754/Muon_Trinity/logs

# ── Environment ───────────────────────────────────────────────────────────────
module purge

# init conda for batch shell, then activate
source /uufs/chpc.utah.edu/common/home/u1520754/miniconda3/etc/profile.d/conda.sh
conda activate jupyter_env

# ── Force CPU (no GPU on this partition) ──────────────────────────────────────
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# ── Sanity check ──────────────────────────────────────────────────────────────
python -c "import torch; print('torch', torch.__version__, 'cuda avail:', torch.cuda.is_available())"

# ── Run ───────────────────────────────────────────────────────────────────────
cd /uufs/chpc.utah.edu/common/home/u1520754/Muon_Trinity/analysis/
python -u train_care_test.py --time_step 10 --agg_method mean

echo "Job finished at $(date)"