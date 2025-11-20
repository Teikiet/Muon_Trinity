import numpy as np
import pickle
import os
import sys

#file = '/scratch/general/nfs1/u6059226/corsika_results/king-mpi/out_sib23d-pId14-En1.00E9-0/OUTPUT_1.00E9_14_0_100/DAT000100-000001'
file = '/uufs/chpc.utah.edu/common/home/u6059226/DAT601285-0001.npy'
#dir_out = '/scratch/general/nfs1/u6059226/corsika_results/king-mpi/out_sib23d-pId14-En1.00E9-0/OUTPUT_1.00E9_14_0_100/'
#if not os.path.exists(dir_out):
#    os.makedirs(dir_out)

data = np.load(file, allow_pickle=True)
print(data)

#muons = []
for pp in data:
    if pp[0] == 13 or pp[0] == -13:
        print(pp[6])
#        muons.append(pp[5])

# Write to pickle file
#with open(os.path.join(dir_out, 'muons_333.pkl'), 'wb') as f:
#    pickle.dump(muons, f)
#print(f'Collected {len(muons)} muons from {file}')


# "pdg": particle ID in PDG scheme (see e.g. https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf)
#  "t": arrival time in ns
#  "x": x position in meter (!)
#  "y": y position in meter (!)
#  "weight": weight factor from thinning
#  "kinetic_energy": kinetic energy GeV
#  "z": production altitude of muons in meter (0 for other particles)
#  "gen": hadron generation counter of muons (0 for other particles)