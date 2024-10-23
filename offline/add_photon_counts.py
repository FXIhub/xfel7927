# doing this first speeds up running
import os
import argparse

if __name__ == '__main__':
    PREFIX = os.environ["EXP_PREFIX"]

    parser = argparse.ArgumentParser(description='add photon counts for a given mask to the events file')
    parser.add_argument('run', type=int, help='Run number')
    parser.add_argument('-n', '--nproc', 
                        help='Number of processes to use',
                        type=int, default=0)
    parser.add_argument('-o', '--out_folder', 
                        help='Path of output folder (default=%s/scratch/events/)'%PREFIX,
                        default=PREFIX+'/scratch/events/')
    parser.add_argument('-m', '--masks', 
                        help=f'pixel masks to apply before calculating litpixels and photon counts located in {PREFIX}scratch/det/. Mask must be 1 for good and 0 for bad, and the data must be in /entry_1/good_pixels)',
                        type=str, nargs='+',
                        default=['r0551_mask.h5'] )
    args = parser.parse_args()
    
    vds_file     = PREFIX+'scratch/vds/r%.4d.cxi' %args.run
    
    # for testing
    #vds_file     = PREFIX+'scratch/saved_hits/r%.4d_hits.cxi' %args.run
    
    events_fname      = args.out_folder + os.path.splitext(os.path.basename(vds_file))[0] + '_events.h5'

import common
from constants import DET_DIST, VDS_DATASET, VDS_MASK_DATASET
import subprocess
import numpy as np
import h5py
import ctypes
import multiprocessing as mp
import utils
from pathlib import Path
import extra_geom
import sys
import time


class Calc_photons():
    def __init__(self, vds_file, mask, nproc=0):
        self.vds_file = vds_file
        if nproc == 0:
            self.nproc = int(subprocess.check_output('nproc').decode().strip())
        else:
            self.nproc = nproc
        print('Using %d processes' % self.nproc)
        
        with h5py.File(vds_file, 'r') as f:
            self.dshape = f[VDS_DATASET].shape
        
        self.mask    = mask
        self.frames  = int(self.dshape[0])
        self.modules = int(self.dshape[1])
        self.pixels  = int(np.prod(self.dshape[2:]))
        self.module_shape = self.dshape[2:]
        self.frame_shape = self.dshape[1:]
        
    def run_frame(self):
        sys.stdout.write('Calculating number of lit pixels for all modules for %d frames\n'%(self.frames))
        sys.stdout.flush()
        
        # radial profile for each frame
        photons   = mp.Array(ctypes.c_ulong, (self.frames,))
        
        jobs = []
        for c in range(self.nproc):
            p = mp.Process(target=self._part_worker, args=(c, photons))
            jobs.append(p)
            p.start()
        
        for i, j in enumerate(jobs):
            sys.stdout.write(f'\n\nwaiting for job {i} to finish... ')
            sys.stdout.flush()
            j.join()
            sys.stdout.write(f'Done\n\n')
            sys.stdout.flush()
        
        self.photons   = np.frombuffer(photons.get_obj(), dtype='u8')
        return self.photons
    
    def _part_worker(self, p, photons):
        np_photons   = np.frombuffer(photons.get_obj(),   dtype='u8')
        
        nframes  = self.frames
        my_start = (nframes // self.nproc) * p
        my_end   = min((nframes // self.nproc) * (p+1), nframes)
        Nevents  = my_end - my_start 
        
        stime = time.time()
        f_vds = h5py.File(self.vds_file, 'r')
        for d in range(my_start, my_end):
            # read vds file (assume photon units)
            mask = self.mask * (f_vds[VDS_MASK_DATASET][d] == 0)
             
            vals = mask * f_vds[VDS_DATASET][d]
            
            photons[d] = np.sum(vals)
            
            etime = time.time()
            if p == 0:
                sys.stdout.write('\r%.4d/%.4d: %.2f Hz' % (d+1, Nevents, (d+1)/(etime-stime)*self.nproc))
                sys.stdout.flush()
        f_vds.close()
        if p == 0:
            sys.stdout.write('\n')
            sys.stdout.flush()


if __name__ == '__main__':
    # check that vds file exists
    d = Path(vds_file).resolve()
    assert(d.is_file())
    
    # check that output directory exists
    d = Path(args.out_folder).resolve()
    assert(d.is_dir())
    
    # check that output events file exists
    d = Path(events_fname).resolve().parent
    assert(d.is_file())

    mask = None
    for fnam in args.masks:
        with h5py.File(fnam, 'r') as f:
            if mask is None :
                mask = f['entry_1/good_pixels'][()]
            else :
                mask *= f['entry_1/good_pixels'][()]
    
    calc_photons = Calc_photons(vds_file, mask, nproc=args.nproc)
    
    photons = calc_photons.run_frame()
    
    note = f"""
    Used {args.masks} for pixel masking in addition to the mask in VDS file
    """
    
    print(f'writing hits powder to {events_fname}')
    with h5py.File(events_fname, 'r+') as f:
        utils.update_h5(f, 'photon_counts', photons, compression=True)
        utils.update_h5(f, 'note', note_hits, compression=False)



