"""
add cellID, pulseID, trainID, detector distance to an events file, eg events/events_r0034.h5
"""

# doing this first speeds up running
import os
import argparse

PREFIX = os.environ["EXP_PREFIX"]

parser = argparse.ArgumentParser(description='Lit pixel calculator')
parser.add_argument('run', type=int, help='Run number')
parser.add_argument('-n', '--nproc', 
                    help='Number of processes to use',
                    type=int, default=0)
parser.add_argument('-o', '--out_folder', 
                    help='Path of output folder (default=%s/scratch/events/)'%PREFIX,
                    default=PREFIX+'/scratch/events/')
parser.add_argument('--out_folder_powder', 
                    help='Path of output folder (default=%s/scratch/powder/)'%PREFIX,
                    default=PREFIX+'/scratch/powder/')
parser.add_argument('-m', '--masks', 
                    help=f'pixel masks to apply before calculating litpixels and photon counts located in {PREFIX}scratch/det/. Mask must be 1 for good and 0 for bad, and the data must be in /entry_1/good_pixels)',
                    type=str, nargs='+',
                    default=['r0065_mask.h5', 'hit_finding_mask.h5'] )
args = parser.parse_args()

args.masks = [f'{PREFIX}scratch/det/{f}' for f in args.masks]

import sys
import time
import glob
import multiprocessing as mp
import ctypes
import subprocess

import h5py
import numpy as np

import utils
from constants import VDS_DATASET, VDS_MASK_DATASET, NMODULES, CHUNK_SIZE, BAD_CELLIDS



class LitPixels():
    def __init__(self, vds_file, mask, nproc=0, chunk_size=CHUNK_SIZE, total_intens=False):
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

        # test
        #self.frames  = 32 * 1000
        
        self.mask    = mask
        self.frames  = int(self.dshape[0])
        self.modules = int(self.dshape[1])
        self.pixels  = int(np.prod(self.dshape[2:]))
        self.module_shape = self.dshape[2:]
        self.frame_shape = self.dshape[1:]

    def run_frame(self):
        sys.stdout.write('Calculating number of lit pixels for all modules for %d frames\n'%(self.frames))
        sys.stdout.flush()
        
        # Litpixels for each module and each frame
        litpix = mp.Array(ctypes.c_ulong, self.frames)
        
        # photon counts for each module and each frame
        intens = mp.Array(ctypes.c_ulong, self.frames)
        
        # powder for each module 
        powder = mp.Array(ctypes.c_ulong, self.modules * self.pixels)
        
        jobs = []
        for c in range(self.nproc):
            p = mp.Process(target=self._part_worker, args=(c, litpix, intens, powder))
            jobs.append(p)
            p.start()
        
        for i, j in enumerate(jobs):
            sys.stdout.write(f'\n\nwaiting for job {i} to finish... ')
            sys.stdout.flush()
            j.join()
            sys.stdout.write(f'Done\n\n')
            sys.stdout.flush()
        
        self.litpix = np.frombuffer(litpix.get_obj(), dtype='u8').reshape(self.frames)
        self.intens = np.frombuffer(intens.get_obj(), dtype='u8').reshape(self.frames)
        self.powder = np.frombuffer(powder.get_obj(), dtype='u8').reshape(self.frame_shape)
        return self.litpix, self.intens, self.powder
    
    def _part_worker(self, p, litpix, intens, powder):
        np_litpix = np.frombuffer(litpix.get_obj(), dtype='u8').reshape(self.frames)
        np_intens = np.frombuffer(intens.get_obj(), dtype='u8').reshape(self.frames)
        np_powder = np.frombuffer(powder.get_obj(), dtype='u8').reshape(self.frame_shape)
        
        nframes = self.frames
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
            mask = f_vds[VDS_MASK_DATASET][pmin:pmax] == 0
            cids = f_vds['entry_1/cellId'][pmin:pmax]
            vals = f_vds[VDS_DATASET][pmin:pmax]
            
            # mask bad cells and pixels
            vals[np.isin(cids, BAD_CELLIDS)] = 0
            vals *= mask
             
            t = np.sum(vals, axis=0, dtype=np_powder.dtype)
            #print(np_powder.dtype, np_powder.shape, vals.dtype, vals.shape, t.dtype, t.shape)
            np_powder           += t
             
            # apply hitfinding mask
            vals *= self.mask
            
            np_litpix[pmin:pmax] = np.sum(vals>0, axis = (1,2,3))
            np_intens[pmin:pmax] = np.sum(vals, axis = (1,2,3))
            
            etime = time.time()
            if p == 0:
                sys.stdout.write('\r%.4d/%.4d: %.2f Hz' % (c+1, num_chunks, (c+1)*self.chunk_size/(etime-stime)*self.nproc))
                sys.stdout.flush()
        f_vds.close()
        if p == 0:
            sys.stdout.write('\n')
            sys.stdout.flush()

vds_file     = PREFIX+'scratch/vds/r%.4d.cxi' %args.run
out_fname    = args.out_folder + os.path.splitext(os.path.basename(vds_file))[0] + '_events.h5'
powder_fname = args.out_folder_powder + os.path.splitext(os.path.basename(vds_file))[0] + '_powder.h5'
modules   = range(NMODULES)

print(f'reading data from {vds_file}')
print(f'output file       {out_fname}')
print(f'powder file       {powder_fname}')

# check that vds file exists
assert(os.path.exists(vds_file))

# check that vds file points to corrected data
result = subprocess.run(['h5ls', '-rv', vds_file], stdout=subprocess.PIPE)
assert('CORR' in result.stdout.decode('utf-8'))

print(f'loading masks {args.masks}')
if len(args.masks) > 0 :
    mask = None
    for fnam in args.masks:
        with h5py.File(fnam) as f:
            if mask is None :
                mask  = f['entry_1/good_pixels'][()]
            else :
                mask *= f['entry_1/good_pixels'][()]
else :
    mask = 1

print('Calculating lit pixels from', vds_file)
l = LitPixels(vds_file, mask, nproc=args.nproc)

print('Running on the following modules:', modules)

# get number of lit pixels, integrated photon counts and powder for each module
litpixels, photon_counts, powder = l.run_frame()

# get train, cell and pulse ID for each frame
print('Copying IDs from VDS file')
sys.stdout.flush()
with h5py.File(vds_file, 'r') as f_vds:
    trainId = f_vds['entry_1/trainId'][:]
    cellId  = f_vds['entry_1/cellId'][:]
    pulseId = f_vds['entry_1/pulseId'][:]

with h5py.File(out_fname, 'a') as f:
    utils.update_h5(f, 'total_intens', np.array(photon_counts), compression=True)
    utils.update_h5(f, 'litpixels', np.array(litpixels), compression=True)
    
    utils.update_h5(f, 'trainId', trainId, compression=True)
    utils.update_h5(f, 'cellId', cellId, compression=True)
    utils.update_h5(f, 'pulseId', pulseId, compression=True)

with h5py.File(powder_fname, 'a') as f:
    utils.update_h5(f, 'powder', np.array(powder), compression=True)
    
print('DONE')
                
