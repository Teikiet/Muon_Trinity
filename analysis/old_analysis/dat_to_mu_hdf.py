from pathlib import Path
import numpy as np
import os
import sys
import pickle
from corsikaio import CorsikaParticleFile

muon_mass = 0.1056583745  # GeV
max_muons = 100000

# Arguments
model = sys.argv[1]
pid = sys.argv[2]
energy = sys.argv[3]
num = sys.argv[4]

# File paths and directory setup
file = f'/scratch/general/nfs1/u6059226/corsika_results/king/out_{model}-pId{pid}-En1.00E{energy}-0/OUTPUT_1.00E{energy}_{pid}_0_{num}/DAT000{num}'
filename = os.path.basename(file)
hdf_dir = f'/scratch/general/nfs1/u6059226/corsika_results/king/hdf/{model}'
data_tag = f'out_{model}-pId{pid}-En1.00E{energy}-0'

print('Using file', file)
print('Max muons per file:', max_muons)
sub_file_counter = 0

# Particle file reading and muon collection
with CorsikaParticleFile(file, parse_blocks=False) as f:
    for event in f:
        print('Start reading particles')
        collected_mu = np.ones(max_muons) * -1  # Initialize with invalid entries
        mu_counter = 0
        for pp in event.particles:
            particle_pid = np.abs(int(pp[0] // 1000))  # Determine particle ID
            if particle_pid in [5, 6]:  # Muon check
                energy = np.sqrt(muon_mass**2 + pp[1]**2 + pp[2]**2 + pp[3]**2)
                collected_mu[mu_counter] = energy
                mu_counter += 1

                # Save batch if max_muons reached
                if mu_counter == max_muons:
                    hdf_file_dir_out = os.path.join(hdf_dir, data_tag + '-' + filename[-3:])
                    os.makedirs(hdf_file_dir_out, exist_ok=True)
                    
                    # Write to pickle file
                    pkl_file_out = os.path.join(hdf_file_dir_out, f'muons_{sub_file_counter}.pkl')
                    with open(pkl_file_out, 'wb') as data_file:
                        pickle.dump(collected_mu, data_file)
                    print(f'Writing subfile {sub_file_counter} to {pkl_file_out}')
                    
                    # Reset for next batch
                    sub_file_counter += 1
                    collected_mu = np.ones(max_muons) * -1
                    mu_counter = 0

        # Save any remaining muons after finishing particle loop
        if mu_counter > 0:  # Only if there are muons collected
            hdf_file_dir_out = os.path.join(hdf_dir, data_tag + '-' + filename[-3:])
            os.makedirs(hdf_file_dir_out, exist_ok=True)
            
            # Write remaining muons to pickle file
            pkl_file_out = os.path.join(hdf_file_dir_out, f'muons_{sub_file_counter}.pkl')
            with open(pkl_file_out, 'wb') as data_file:
                pickle.dump(collected_mu[:mu_counter], data_file)  # Only write collected muons
            print(f'Writing final subfile {sub_file_counter} to {pkl_file_out}')
