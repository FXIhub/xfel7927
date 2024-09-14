import argparse
import os
import sys
import time

import writeemc
import h5py
import numpy as np

from constants import PREFIX

parser = argparse.ArgumentParser(description='Save hits to emc file')
parser.add_argument('run', help='Run number', type=int)
args = parser.parse_args()

cxi_file = PREFIX + '/scratch/saved_hits/r%.4d_hits.cxi'%args.run

# Merge modules
wemc_all = writeemc.EMCWriter(PREFIX+'/scratch/emc/r%.4d.emc' % args.run, 1024**2, hdf5=False)
stime = time.time()

with h5py.File(cxi_file) as f:
    photons = f['/entry_1/data_1/data']
    
    for i in range(photons.shape[0]):
        wemc_all.write_frame(photons[i].ravel())
        
        if (i+1) % 10 == 0:
            sys.stderr.write('\rWritten frame %d/%d (%.3f Hz)' % (i+1, photons.shape[0], (i+1)/(time.time()-stime)))
            sys.stderr.flush()
sys.stderr.write('\n')
sys.stderr.flush()

wemc_all.finish_write()

