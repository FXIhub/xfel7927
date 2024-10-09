import argparse
from constants import PREFIX, DET_DIST, FRAME_SHAPE, VDS_DATASET, VDS_MASK_DATASET
import common

parser = argparse.ArgumentParser(description='Calculate average background from non-hits and add to a cxi file')
parser.add_argument('run', type=int, help='Run number')

args = parser.parse_args()
    
args.output_file   = PREFIX+'scratch/saved_hits/r%.4d_hits.cxi' %args.run
args.vds_file      = PREFIX+'scratch/vds/r%.4d.cxi' %args.run
args.events_file   = PREFIX+'scratch/events/r%.4d_events.h5'%args.run

import numpy as np
import h5py
#import extra_data
from tqdm import tqdm
import sys
import time
import utils

import multiprocessing as mp

"""
h5ls r0035_events.h5
/cellId                  Dataset {1055808, 16}
/is_hit                  Dataset {1055808}
/is_miss                 Dataset {1055808}
/litpixels               Dataset {1055808, 16}
/pulseId                 Dataset {1055808}
/pulse_energy            Dataset {1055808}
/total_intens            Dataset {1055808, 16}
/trainId                 Dataset {1055808}
/wavelength              Dataset {1055808}
"""

print(f'loading miss selection from {args.events_file}')
with h5py.File(args.events_file) as f:
    # indices for definite miss
    #m = f['is_miss'][()]
    m = ~f['is_hit'][()]
    indices = np.where(m)[0]
    
    if len(f['cellId'].shape) == 2 :
        cellId_lit    = f['cellId'][:, 0]
    elif len(f['cellId'].shape) == 1 :
        cellId_lit    = f['cellId'][:]
    else :
        raise ValueError(f'unknown cellId shape in {args.events_file}')

    trainId_lit   = f['trainId'][()]
    pulseId_lit   = f['pulseId'][()]

Nevents = len(indices)

print(f'found {Nevents} events labeled as miss')
sys.stdout.flush()

#size = mp.cpu_count()
size = 16

# split frames over ranks
events_rank = np.linspace(0, Nevents, size+1).astype(int)

# load event ids for cross-checking
with h5py.File(args.vds_file) as g:
    cellId_vds = g['entry_1/cellId'][:, 0]
    trainId_vds = g['entry_1/trainId'][()]
    pulseId_vds = g['entry_1/pulseId'][()]

def worker(rank, b):
    my_indices = indices[events_rank[rank]: events_rank[rank+1]] 
    
    print(f'rank {rank} is processing indices {events_rank[rank]} to {events_rank[rank+1]}')
    sys.stdout.flush()

    if rank == 0 :
        it = tqdm(range(len(my_indices)), desc = f'Processing data from {args.vds_file}')
    else :
        it = range(len(my_indices))

    background = np.zeros(FRAME_SHAPE, dtype=float)
    overlap    = np.zeros(FRAME_SHAPE, dtype=int)
    
    with h5py.File(args.vds_file) as g:
        data = g[VDS_DATASET]
        mask = g[VDS_MASK_DATASET]
        
        for i in it:
            index = my_indices[i]
            
            frame = np.squeeze(data[index])
            m     = np.squeeze(mask[index] == 0)
            
            background += frame * m
            overlap    += m
            
            # make sure we have the right event
            assert(cellId_vds[index] == cellId_lit[index])
            assert(trainId_vds[index] == trainId_lit[index])
            assert(pulseId_vds[index] == pulseId_lit[index])

    b.put([background, overlap])


b = mp.Queue()
jobs = [mp.Process(target=worker, args=(m, b)) for m in range(size)]
[j.start() for j in jobs]

background_g = np.zeros(FRAME_SHAPE, dtype=float)
overlap_g    = np.zeros(FRAME_SHAPE, dtype=int)
for r in range(size):
    background, overlap = b.get()
    background_g += background
    overlap_g    += overlap
[j.join() for j in jobs]

out = background_g / np.clip(overlap_g, 1, None)

with h5py.File(args.output_file, 'a') as f:
    key = 'entry_1/instrument_1/detector_1/background'
    utils.update_h5(f, key, out.astype(np.float32), compression = True)

print('Done')
sys.stdout.flush()


