#!/bin/bash
#SBATCH --job-name=submitter
#SBATCH --time=70:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=kingspeak-guest
#SBATCH --account=owner-guest
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1520754/logs/submitter_%j.out
#SBATCH --error=/uufs/chpc.utah.edu/common/home/u1520754/logs/submitter_%j.err


bash $HOME/Muon_Trinity/cluster_corsika8/base_muon_event_generation/submit_multiple_continuous_tree.sh --rerun-event
