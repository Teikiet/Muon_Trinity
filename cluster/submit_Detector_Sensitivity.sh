#!/bin/bash

# Array of production heights to test

# Other parameters
#6, 1500:  5, 15000; 4, 15000; 3, 15000
gen_radi="150"
particle=6
N=5004
name_tag="1E3_1E7"  # Change this to your desired name tag
# Create logs directory if it doesn't exist
mkdir -p /uufs/chpc.utah.edu/common/home/u1520754/test_jupyter/submits
mkdir -p /uufs/chpc.utah.edu/common/home/u1520754/test_jupyter/outputs

job_name="detector_pid${particle}_N${N}_${name_tag}"
mkdir -p /uufs/chpc.utah.edu/common/home/u1520754/test_jupyter/submits/pid${particle}_N${N}
# Create individual SLURM script for each height
cat << EOF > "/uufs/chpc.utah.edu/common/home/u1520754/test_jupyter/submits/pid${particle}_N${N}/job_${job_name}.slurm"
#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --account=owner-guest
#SBATCH --partition=kingspeak-shared-guest
#SBATCH --ntasks=1
#SBATCH --mem=16GB
#SBATCH --time=2:00:00
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1520754/test_jupyter/submits/pid${particle}_N${N}/${job_name}.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/u1520754/test_jupyter/submits/pid${particle}_N${N}/${job_name}.err


# Activate your virtual environment if you have one
source /uufs/chpc.utah.edu/common/home/u1520754/miniconda3/etc/profile.d/conda.sh
conda activate jupyter_env

# Run the simulation
cd /uufs/chpc.utah.edu/common/home/u1520754/cluster
python Detector_Sensitivity_MonteCarlos.py \\
    --particle ${particle} \\
    --gen_radi ${gen_radi} \\
    --N ${N} \\
    --noise ${noise}\\
    --job_id \${SLURM_JOB_ID} \\
    --name_tag ${name_tag}

echo "Job completed "
EOF

# Submit the job
sbatch "/uufs/chpc.utah.edu/common/home/u1520754/test_jupyter/submits/pid${particle}_N${N}/job_${job_name}.slurm"
echo "All jobs submitted!"
