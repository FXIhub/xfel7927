
# doing this first speeds up running
import os
import argparse

if __name__ == '__main__':
    PREFIX = os.environ["EXP_PREFIX"]

    parser = argparse.ArgumentParser(description='calculate radial profiles for every frame')
    parser.add_argument('run', type=int, help='Run number')
    parser.add_argument('-n', '--nproc', 
                        help='Number of processes to use',
                        type=int, default=0)
    parser.add_argument('-o', '--out_folder', 
                        help='Path of output folder (default=%s/scratch/radial_profiles/)'%PREFIX,
                        default=PREFIX+'/scratch/radial_profiles/')
    parser.add_argument('-m', '--masks', 
                        help=f'pixel masks to apply before calculating litpixels and photon counts located in {PREFIX}scratch/det/. Mask must be 1 for good and 0 for bad, and the data must be in /entry_1/good_pixels)',
                        type=str, nargs='+',
                        default=['r0551_mask.h5'] )
    args = parser.parse_args()

    vds_file     = PREFIX+'scratch/vds/r%.4d.cxi' %args.run
    
    # for testing
    #vds_file     = PREFIX+'scratch/saved_hits/r%.4d_hits.cxi' %args.run

    out_fname    = args.out_folder + os.path.splitext(os.path.basename(vds_file))[0] + '_radial_profiles.h5'
    for i in range(len(args.masks)):
        args.masks[i] = f'{PREFIX}/scratch/det/{args.masks[i]}'


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


class Radial_average():
    def __init__(self, xyz, mask, radial_bin_size):
        self.r = (xyz[0]**2 + xyz[1]**2)**0.5
        self.radial_bin_size = radial_bin_size
        
        self.mask = mask.copy()
        
        # integer label for radial bin
        self.rl = np.rint(self.r[self.mask] / radial_bin_size).astype(int).ravel()
        
        # just to speed things up a little
        self.rl = np.ascontiguousarray(self.rl)
        
        # number of pixels contributing to each bin 
        self.rcounts = np.bincount(self.rl.ravel())
        
        # for reference
        self.rs = np.arange(np.max(self.rl)+1) * self.radial_bin_size
        
    def make_rad_av(self, ar):
        rsum = np.bincount(self.rl, ar[self.mask])
        
        # normalise, might not want to do this
        rsum /= np.clip(self.rcounts, 1, None)
        
        # to broadcast back to image
        #im = np.zeros(self.mask.shape, dtype = float)
        #im[self.mask] = rsum[self.rl]
        return rsum


class Calc_rad():
    def __init__(self, vds_file, mask, xyz, radial_bin_size = 200e-6, nproc=0):
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
        
        # initialise radial profile calculation 
        self.ra = Radial_average(xyz, mask, radial_bin_size)
        self.radial_pixels = int(self.ra.rs.shape[0])
    
    def run_frame(self):
        sys.stdout.write('Calculating number of lit pixels for all modules for %d frames\n'%(self.frames))
        sys.stdout.flush()
        
        # radial profile for each frame
        radial_profiles = mp.Array(ctypes.c_float, self.frames * self.radial_pixels)
        
        jobs = []
        for c in range(self.nproc):
            p = mp.Process(target=self._part_worker, args=(c, radial_profiles))
            jobs.append(p)
            p.start()
        
        for i, j in enumerate(jobs):
            sys.stdout.write(f'\n\nwaiting for job {i} to finish... ')
            sys.stdout.flush()
            j.join()
            sys.stdout.write(f'Done\n\n')
            sys.stdout.flush()
        
        self.radial_profiles = np.frombuffer(radial_profiles.get_obj(), dtype=np.float32).reshape((self.frames, self.radial_pixels))
        return self.radial_profiles
    
    def _part_worker(self, p, radial_profiles):
        np_radial_profiles = np.frombuffer(radial_profiles.get_obj(), dtype=np.float32).reshape((self.frames, self.radial_pixels))
        
        nframes  = self.frames
        my_start = (nframes // self.nproc) * p
        my_end   = min((nframes // self.nproc) * (p+1), nframes)
        Nevents  = my_end - my_start 

        stime = time.time()
        f_vds = h5py.File(self.vds_file, 'r')
        for d in range(my_start, my_end):
            # read vds file (assume photon units)
            mask = self.mask * (f_vds[VDS_MASK_DATASET][d] == 0)
             
            # for testing
            #mask = self.mask

            cids = f_vds['entry_1/cellId'][d]
            vals = f_vds[VDS_DATASET][d]
            
            np_radial_profiles[d] = self.ra.make_rad_av(vals)
            
            etime = time.time()
            if p == 0:
                sys.stdout.write('\r%.4d/%.4d: %.2f Hz' % (d+1, Nevents, (d+1)/(etime-stime)*self.nproc))
                sys.stdout.flush()
        f_vds.close()
        if p == 0:
            sys.stdout.write('\n')
            sys.stdout.flush()


if __name__ == '__main__':
    # get geom object
    args.geom_file = common.get_geom(args.run)
    geom           = extra_geom.AGIPD_1MGeometry.from_crystfel_geom(args.geom_file)
    xyz            = np.transpose(geom.get_pixel_positions(), (3, 0, 1, 2))
    
    # check that vds file exists
    d = Path(vds_file).resolve()
    assert(d.is_file())
    
    # check that output directory exists
    d = Path(out_fname).resolve().parent
    assert(d.is_dir())
    
    mask = np.ones(xyz.shape[1:], dtype = bool)
    for fnam in args.masks:
        with h5py.File(fnam, 'r') as f:
            mask *= f['entry_1/good_pixels'][()]
    
    calc_rad = Calc_rad(vds_file, mask, xyz, radial_bin_size = 200e-6, nproc=args.nproc)
    
    radial_profiles = calc_rad.run_frame()

    note = """
    To broadcast the radial averages back to the frame geometry:
        frame = np.zeros(mask.shape, dtype = float)
        frame[mask] = radial_profiles[0][radial_bin_labels]
    """
    
    print(f'writing radial profiles to {out_fname}')
    with h5py.File(out_fname, 'a') as f:
        utils.update_h5(f, 'data', radial_profiles, compression=True)
        utils.update_h5(f, 'radial_bin_labels',   calc_rad.ra.rl, compression=True)
        utils.update_h5(f, 'radial_bin_counts',   calc_rad.ra.rcounts, compression=True)
        utils.update_h5(f, 'radial_pixel_values', calc_rad.ra.rs, compression=True)
        utils.update_h5(f, 'mask', calc_rad.ra.mask, compression=True)
        utils.update_h5(f, 'note', note, compression=False)
