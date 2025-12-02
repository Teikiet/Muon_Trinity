#!/bin/bash
#SBATCH --account=owner-guest
#SBATCH --partition=kingspeak-shared-guest
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --job-name=corsika_shower
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1520754/corsika8_output/logs/corsika_shower_%j.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/u1520754/corsika8_output/logs/corsika_shower_%j.err

# Print job info
echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Started at: $(date)"
echo "Working directory: $(pwd)"
echo "=================================================="

# Load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate corsika

# Add CORSIKA to PATH
export PATH=$HOME/corsika-install/bin:$PATH

# Create output directories
OUTPUT_DIR=$HOME/corsika8_output
LOG_DIR=$OUTPUT_DIR/logs
mkdir -p $OUTPUT_DIR
mkdir -p $LOG_DIR

# Navigate to output directory (important!)
cd $OUTPUT_DIR
rm -rf my_shower

echo "Output directory: $OUTPUT_DIR"
echo "Log directory: $LOG_DIR"
echo ""

# Run CORSIKA air shower simulation
# Parameters:
#   --pdg 2212   : Proton primary particle
#   -E 1e5       : Energy in GeV (100 TeV)
#   -f my_shower : Output filename prefix

echo "Running CORSIKA simulation..."
echo "Primary: Proton"
echo "Energy: 1e5 GeV (100 TeV)"
echo "Output prefix: my_shower"
echo ""

# Run and capture exit status
c8_air_shower --pdg 2212 -E 1e5 -f my_shower
EXIT_STATUS=$?

# Check if simulation completed
if [ $EXIT_STATUS -eq 0 ]; then
    echo ""
    echo "=================================================="
    echo "Simulation completed successfully!"
    echo "Output files in: $OUTPUT_DIR"
    echo ""
    ls -lh my_shower* 2>/dev/null || echo "No output files found!"
    echo ""
    echo "Total disk usage:"
    du -sh $OUTPUT_DIR
    echo "=================================================="
else
    echo ""
    echo "=================================================="
    echo "ERROR: Simulation failed with exit code $EXIT_STATUS"
    echo "Check error log at: $LOG_DIR/corsika_shower_${SLURM_JOB_ID}.err"
    echo "=================================================="
    exit $EXIT_STATUS
fi

echo "Finished at: $(date)"
