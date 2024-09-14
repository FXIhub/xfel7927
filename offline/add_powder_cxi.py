import argparse
from constants import PREFIX, DET_DIST, FRAME_SHAPE, VDS_DATASET, VDS_MASK_DATASET
import common

parser = argparse.ArgumentParser(description='Calculate sum of hits and write to cxi file')
parser.add_argument('run', type=int, help='Run number')

args = parser.parse_args()
    
args.output_file   = PREFIX+'scratch/saved_hits/r%.4d_hits.cxi' %args.run

import numpy as np
import h5py
#import extra_data
from tqdm import tqdm
import sys
import time
import utils

import multiprocessing as mp

print(f'Calculating powder pattern of hits from {args.output_file}')

with h5py.File(args.output_file) as f:
    indices = list(range(f['entry_1/data_1/data'].shape[0]))
    
Nevents = len(indices)

size = 16

# split frames over ranks
events_rank = np.linspace(0, Nevents, size+1).astype(int)

def worker(rank, b):
    my_indices = indices[events_rank[rank]: events_rank[rank+1]] 
    
    print(f'rank {rank} is processing indices {events_rank[rank]} to {events_rank[rank+1]}')
    sys.stdout.flush()
    
    if rank == 0 :
        it = tqdm(range(len(my_indices)), desc = f'Processing data from {args.output_file}')
    else :
        it = range(len(my_indices))

    powder = np.zeros(FRAME_SHAPE, dtype=float)
    overlap    = np.zeros(FRAME_SHAPE, dtype=int)
    
    with h5py.File(args.output_file) as g:
        data = g[VDS_DATASET]
        mask = g[VDS_MASK_DATASET]
        
        for i in it:
            index = my_indices[i]
            
            frame = np.squeeze(data[index])
            m     = np.squeeze(mask[index] == 0)
            
            powder     += frame * m
            overlap    += m
            
    b.put([powder, overlap])


b = mp.Queue()
jobs = [mp.Process(target=worker, args=(m, b)) for m in range(size)]
[j.start() for j in jobs]

powder_g = np.zeros(FRAME_SHAPE, dtype=float)
overlap_g    = np.zeros(FRAME_SHAPE, dtype=int)
for r in range(size):
    powder, overlap = b.get()
    powder_g += powder
    overlap_g    += overlap
[j.join() for j in jobs]

out = powder_g / np.clip(overlap_g, 1, None)

with h5py.File(args.output_file, 'a') as f:
    key = 'entry_1/instrument_1/detector_1/powder'
    utils.update_h5(f, key, out.astype(np.float32), compression = True)

print('Done')
sys.stdout.flush()
