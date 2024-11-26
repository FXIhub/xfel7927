import argparse
from constants import PREFIX, DET_DIST, FRAME_SHAPE, VDS_DATASET, VDS_MASK_DATASET
import common

parser = argparse.ArgumentParser(description='Calculate average background from non-hits and add to a cxi file')
parser.add_argument('run', type=int, help='Run number')

args = parser.parse_args()
    
args.output_file   = PREFIX+'scratch/saved_hits/r%.4d_hits.cxi' %args.run
args.vds_file      = PREFIX+'scratch/vds/r%.4d.cxi' %args.run
args.events_file   = PREFIX+'scratch/events/r%.4d_events.h5'%args.run
args.back_file     = PREFIX+'scratch/powder/r%.4d_powder_is_hit_False_per_pixel.h5'%args.run

import numpy as np
import h5py
#import extra_data
from tqdm import tqdm
import sys
import time
import utils


"""
h5ls -r r0600_events.h5
/                        Group
/cellId                  Dataset {1052649}
/hit_median_train        Dataset {1052649}
/hit_sigma               Dataset {1052649}
/hit_std_train           Dataset {1052649}
/is_hit                  Dataset {1052649}
/is_miss                 Dataset {1052649}
/litpixels               Dataset {1052649}
/litpixels_mask          Dataset {1052649}
/modules                 Dataset {16}
/pulse_energy            Dataset {1052649}
/total_intens            Dataset {1052649}
/total_intens_mask       Dataset {1052649}
/trainId                 Dataset {1052649}
/wavelength              Dataset {1052649}


h5ls -r r0600_powder_is_hit_False_per_pixel.h5
/                        Group
/data                    Dataset {16, 512, 128}
"""

# calculate background adjustment factor per train
#
# a_t = sum_(d misses in train t) K_di / (<B>_i x number of misses in train t)
#
# b_d = e_d x a_d / <e> 
# where e_d is the pulse energy and a_(t_d) is the adjustment factor for frame d

# minimum reliable energy reading
emin = 1e-3

# mean pulse energy 
with h5py.File(args.events_file) as f:
    e      = f['pulse_energy'][()]
    
    trainIds = f['/trainId'][()]
    misses   = ~f['/is_hit'][()]
    photon_counts = f['/total_intens'][()]

unique_trainIds = np.unique(trainIds)

# mean background 
with h5py.File(args.back_file) as f:
    back = f['data'][()]

# mean background signal per shot
back_counts = np.sum(back)

# calculate adjustment factor
a_d = -np.ones(trainIds.shape[0], dtype = float)

for train in unique_trainIds:
    m = np.where(train == trainIds)[0]

    # misses in train
    n = m[misses[m]]
    
    if len(n) > 0 :
        # sum K
        Ksum = np.sum(photon_counts[n])
        
        a_t = Ksum / (back_counts * len(n))
        a_d[m] = a_t

# b = e_d a_d / <e>
b = -np.ones(a_d.shape[0], dtype = float)

# don't rely on energy readings above 1e-3
m = e > emin
e[m]  = e[m] / np.mean(e[m])
e[~m] = 1

b = e * a_d 

# write to events file
with h5py.File(args.events_file, 'r+') as f:
    key = 'background_weighting'
    utils.update_h5(f, key, b.astype(np.float32), compression = True)

# write to cxi file
with h5py.File(args.output_file, 'a') as f:
    key = 'entry_1/instrument_1/detector_1/background'
    utils.update_h5(f, key, back.astype(np.float32), compression = True)

    vds_index = f['/entry_1/vds_index'][()]
    
    key = '/entry_1/background_weighting'
    utils.update_h5(f, key, b[vds_index].astype(np.float32), compression = True)


print('add_background Done')
sys.stdout.flush()


