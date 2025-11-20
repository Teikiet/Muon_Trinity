#!/bin/bash
#SBATCH --account=owner-guest
#SBATCH --partition=kingspeak-shared-guest
#SBATCH --job-name=corsika
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1520754/example_out_88/example.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/u1520754/example_out_88/example.err
#SBATCH --mem=6G           # Request less memory
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mail-user=u1520754@utah.edu

module load gcc
cd /uufs/chpc.utah.edu/common/home/u1520754/corsika-78010/run
ulimit -s unlimited
rm /uufs/chpc.utah.edu/common/home/u1520754/example_out_88/*
mkdir -p /uufs/chpc.utah.edu/common/home/u1520754/example_out_88

EXEC=./corsika78010Linux_SIBYLL_urqmd
INPUT=/uufs/chpc.utah.edu/common/home/u1520754/cluster/example.inp

if [ ! -f "$EXEC" ]; then
  echo "ERROR: Executable not found: $EXEC"
  exit 1
fi

if [ ! -f "$INPUT" ]; then
  echo "ERROR: Input file not found: $INPUT"
  exit 1
fi

$EXEC < $INPUT > /uufs/chpc.utah.edu/common/home/u1520754/example_out_88/example.log