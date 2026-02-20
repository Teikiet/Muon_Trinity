#!/bin/bash
#SBATCH --account=owner-guest
#SBATCH --partition=kingspeak-guest
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --job-name=corsika_shower
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1520754/corsika8_output_4/logs/corsika_shower_muon.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/u1520754/corsika8_output_4/logs/corsika_shower_muon.err

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
OUTPUT_DIR=$HOME/corsika8_output_4
LOG_DIR=$OUTPUT_DIR/logs
mkdir -p $OUTPUT_DIR
mkdir -p $LOG_DIR

# Navigate to output directory (important!)
cd $OUTPUT_DIR
#name of output files
OUTPUT_FILE="muon_test"
rm -rf ${OUTPUT_FILE}*

echo "Output directory: $OUTPUT_DIR"
echo "Log directory: $LOG_DIR"
echo ""

# Run CORSIKA air shower simulation

echo "Running CORSIKA simulation..."
echo "Output prefix: $OUTPUT_FILE"
echo ""
#13(mu-) -13(mu+) 15(tau-), -15(tau+) 2212 is proton seed with muon shower: 2, 5, 6 at 20km injection height
# Run and capture exit status
muon_cherenkov --help
muon_cherenkov --pdg 15 -E 1e6 \
 --injection-height 3000 \
 --emcut 0.02 --mucut 0.43 --hadcut 0.5 \
 --observation-level 2944 --zenith 89 --azimuth 270 \
 --cherenkov-projection telescope --telescope-azimuth 270 --telescope-zenith 89 \
 --telescope-x 0 --telescope-y 0 --telescope-z 0 --telescope-radius 5 \
 --seed 5 \
 --hadronModel SIBYLL-2.3d -f $OUTPUT_FILE \
 --cherenkov-output-format eventio \
 --atmosphere-absorption-file /uufs/chpc.utah.edu/common/home/u1520754/corsika/modules/data/CHERENKOV/atmosphere/abstable_MODTRAN_new.dat \
 #--enable-curvature \
 #--disable-dispersion

EXIT_STATUS=$?

# Check if simulation completed
if [ $EXIT_STATUS -eq 0 ]; then
    echo ""
    echo "=================================================="
    echo "Simulation completed successfully!"
    echo "Output files in: $OUTPUT_DIR"
    echo ""
    ls -lh ${OUTPUT_FILE}* 2>/dev/null || echo "No output files found!"
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
