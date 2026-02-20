#!/bin/bash
#SBATCH --account=owner-guest
#SBATCH --partition=kingspeak-guest
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1520754/slurm_log/single_sib23d_slurm-%j.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/u1520754/slurm_log/single_sib23d_slurm-%j.err
#SBATCH --time=4:00:00 
#SBATCH --ntasks=1
#SBATCH --mem=16G


module load gcc
input=${1}
output=${2}

# Print the variables to confirm they are set correctly
echo "Input file: $input"
echo "Output file: $output"

# Navigate to the CORSIKA run directory
cd /uufs/chpc.utah.edu/common/home/u1520754/corsika-78010/run/

# Run the CORSIKA executable with the input file, redirecting the output to the log file
./corsika78010Linux_SIBYLL_urqmd < $input > $output
#./corsika78010Linux_SIBYLL_urqmd_inclined < $input > $output
#sacct -j $SLURM_JOB_ID --format=JobID,JobName,Elapsed,TotalCPU,MaxDiskWrite,MaxRSS > /scratch/general/nfs1/u6059226/corsika_out/sib23d/slurm_cpu_hours_sibyll_${SLURM_JOB_ID}.txt
