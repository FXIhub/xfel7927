"""
add cellID, pulseID, trainID, detector distance to an events file, eg events/events_r0034.h5
"""

# doing this first speeds up running
import os
import argparse

PREFIX = os.environ("EXP_PREFIX")

parser = argparse.ArgumentParser(description='Lit pixel calculator')
parser.add_argument('run', type=int, help='Run number')
parser.add_argument('-d', '--dark_run', type=int, help='Dark run number', default=-1)
parser.add_argument('-n', '--nproc', 
                    help='Number of processes to use',
                    type=int, default=0)
parser.add_argument('-o', '--out_folder', 
                    help='Path of output folder (default=%s/events/)'%PREFIX,
                    default=PREFIX+'/events/')
args = parser.parse_args()

import sys
import time
import glob
import multiprocessing as mp
import ctypes
import subprocess

import h5py
import numpy as np

import utils
from constants import VDS_DATASET, NMODULES, CHUNK_SIZE

print(f'reading data from {vds_file}')
print(f'output file       {out_fname}')


class LitPixels():
    def __init__(self, vds_file, nproc=0, chunk_size=CHUNK_SIZE, total_intens=False):
        self.vds_file = vds_file
        self.chunk_size = chunk_size # Needs to be multiple of 32 for raw data
        if self.chunk_size % 32 != 0:
            print('WARNING: Performance is best with a multiple of 32 chunk_size')
        if nproc == 0:
            self.nproc = int(subprocess.check_output('nproc').decode().strip())
        else:
            self.nproc = nproc
        print('Using %d processes' % self.nproc)

        with h5py.File(vds_file, 'r') as f:
            self.dshape = f[VDS_DATASET].shape

    def run_module(self, module):
        sys.stdout.write('Calculating number of lit pixels in module %d for %d frames\n'%(module, self.dshape[0]))
        sys.stdout.flush()
        
        # Litpixels for each module and each frame
        litpix = mp.Array(ctypes.c_ulong, self.dshape[0])
        
        # photon counts for each module and each frame
        intens = mp.Array(ctypes.c_ulong, self.dshape[0])
        
        # powder for each module 
        powder = mp.Array(ctypes.c_ulong, self.dshape[1:])
        
        jobs = []
        for c in range(self.nproc):
            p = mp.Process(target=self._part_worker, args=(c, module, litpix, intens, powder))
            jobs.append(p)
            p.start()
        
        for j in jobs:
            j.join()
        
        self.litpix = np.frombuffer(litpix.get_obj(), dtype='u8')
        self.intens = np.frombuffer(intens.get_obj(), dtype='u8')
        self.powder = np.frombuffer(intens.get_obj(), dtype='u8')
        return self.litpix, self.intens, self.powder
    
    def _part_worker(self, p, m, litpix, intens, self.powder):
        np_litpix = np.frombuffer(litpix.get_obj(), dtype='u8')
        np_intens = np.frombuffer(intens.get_obj(), dtype='u8')
        np_powder = np.frombuffer(powder.get_obj(), dtype='u8')
        
        nframes = self.dshape[0]
        my_start = (nframes // self.nproc) * p
        my_end = min((nframes // self.nproc) * (p+1), nframes)
        num_chunks = int(np.ceil((my_end-my_start)/self.chunk_size))
        #num_chunks = 4
        if p == 0:
            print('Doing %d chunks of %d frames each' % (num_chunks, self.chunk_size))
            sys.stdout.flush()

        stime = time.time()
        f_vds = h5py.File(self.vds_file, 'r')
        for c in range(num_chunks):
            pmin = my_start + c*self.chunk_size
            pmax = min(my_start + (c+1) * self.chunk_size, my_end)

            # read vds file (assume photon units)
            mask = f_vds[VDS_MASK_DATASET][pmin:pmax, m] == 0
            cids = f_vds['entry_1/cellId'][pmin:pmax, m]
            vals = f_vds[VDS_DATASET][pmin:pmax, m]
            
            # mask bad cells and pixels
            vals[np.isin(cids, BAD_CELLIDS)] = 0
            vals *= mask
            
            np_litpix[pmin:pmax] = np.sum(vals>0, axis = (1,2))
            np_intens[pmin:pmax] = np.sum(vals, axis = (1,2))
            np_powder           += vals
            
            etime = time.time()
            if p == 0:
                sys.stdout.write('\r%.4d/%.4d: %.2f Hz' % (c+1, num_chunks, (c+1)*self.chunk_size/(etime-stime)*self.nproc))
                sys.stdout.flush()
        f_vds.close()
        if p == 0:
            sys.stdout.write('\n')
            sys.stdout.flush()

vds_file  = PREFIX+'vds/r%.4d.cxi' %args.run
out_fname = args.out_folder + os.path.splitext(os.path.basename(vds_file))[0] + '_events.h5'
modules   = range(NMODULES)

print('Calculating lit pixels from', vds_file)
l = LitPixels(vds_file, nproc=args.nproc)

print('Running on the following modules:', modules)

# get number of lit pixels, integrated photon counts and powder for each module
litpixels = []
photon_counts = []
powder = []
for module in modules:
    l, c, p = l.run_module(module)
    litpixels.append(l.copy())
    photon_counts.append(c.copy())
    powder.append(p.copy())

# get train, cell and pulse ID for each frame
print('Copying IDs from VDS file')
sys.stdout.flush()
with h5py.File(vds_fname, 'r') as f_vds:
    trainId = f_vds['entry_1/trainId'][:]
    cellId  = f_vds['entry_1/cellId'][:]
    pulseId = f_vds['entry_1/pulseId'][:]

with h5py.File(out_fname, 'a') as f:
    utils.update_h5(f, 'total_intens', np.array(photon_counts))
    utils.update_h5(f, 'litpixels', np.array(litpixels))
    utils.update_h5(f, 'powder', np.array(powder))
    utils.update_h5(f, 'trainId', trainId)
    utils.update_h5(f, 'cellId', cellId)
    utils.update_h5(f, 'pulseId', pulseId)
    
print('DONE')
                
