# doing this first speeds up running
import os
import argparse

if __name__ == '__main__':
    PREFIX = os.environ["EXP_PREFIX"]

    parser = argparse.ArgumentParser(description='calculate powder profiles for all hits and non hits')
    parser.add_argument('run', type=int, help='Run number')
    parser.add_argument('-n', '--nproc', 
                        help='Number of processes to use',
                        type=int, default=0)
    parser.add_argument('-o', '--out_folder', 
                        help='Path of output folder (default=%s/scratch/powder/)'%PREFIX,
                        default=PREFIX+'/scratch/powder/')
    args = parser.parse_args()
    
    vds_file     = PREFIX+'scratch/vds/r%.4d.cxi' %args.run
    
    # for testing
    #vds_file     = PREFIX+'scratch/saved_hits/r%.4d_hits.cxi' %args.run
    
    events_fname      = args.out_folder + os.path.splitext(os.path.basename(vds_file))[0] + '_events.h5'
    out_fname_hits    = args.out_folder + os.path.splitext(os.path.basename(vds_file))[0] + '_powder_hits.h5'
    out_fname_misses  = args.out_folder + os.path.splitext(os.path.basename(vds_file))[0] + '_powder_misses.h5'

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


class Calc_powder():
    def __init__(self, vds_file, is_hit, is_miss, nproc=0):
        self.vds_file = vds_file
        if nproc == 0:
            self.nproc = int(subprocess.check_output('nproc').decode().strip())
        else:
            self.nproc = nproc
        print('Using %d processes' % self.nproc)
        
        with h5py.File(vds_file, 'r') as f:
            self.dshape = f[VDS_DATASET].shape
        
        self.is_hit  = self.is_hit
        self.is_miss = self.is_miss
        
        self.frames  = int(self.dshape[0])
        self.modules = int(self.dshape[1])
        self.pixels  = int(np.prod(self.dshape[2:]))
        self.module_shape = self.dshape[2:]
        self.frame_shape = self.dshape[1:]
        
    def run_frame(self):
        sys.stdout.write('Calculating number of lit pixels for all modules for %d frames\n'%(self.frames))
        sys.stdout.flush()
        
        # radial profile for each frame
        overlap_hits   = mp.Array(ctypes.c_ulong, self.pixels)
        overlap_misses = mp.Array(ctypes.c_ulong, self.pixels)
        powder_hits    = mp.Array(ctypes.c_ulong, self.pixels)
        powder_misses  = mp.Array(ctypes.c_ulong, self.pixels)
        
        jobs = []
        for c in range(self.nproc):
            p = mp.Process(target=self._part_worker, args=(c, powder_hits, powder_misses))
            jobs.append(p)
            p.start()
        
        for i, j in enumerate(jobs):
            sys.stdout.write(f'\n\nwaiting for job {i} to finish... ')
            sys.stdout.flush()
            j.join()
            sys.stdout.write(f'Done\n\n')
            sys.stdout.flush()
        
        self.overlap_hits   = np.frombuffer(overlap_hits.get_obj(),   dtype='u8').reshape(self.frame_shape)
        self.overlap_misses = np.frombuffer(overlap_misses.get_obj(), dtype='u8').reshape(self.frame_shape)
        self.powder_hits   = np.frombuffer(powder_hits.get_obj(),   dtype='u8').reshape(self.frame_shape)
        self.powder_misses = np.frombuffer(powder_misses.get_obj(), dtype='u8').reshape(self.frame_shape)
        return self.powder_hits, self.powder_misses, self.overlap_hits, self.overlap_misses
    
    def _part_worker(self, p, powder_hits, powder_misses):
        np_powder_hits   = np.frombuffer(powder_hits.get_obj(),   dtype='u8').reshape(self.frame_shape)
        np_powder_misses = np.frombuffer(powder_misses.get_obj(), dtype='u8').reshape(self.frame_shape)
        np_overlap_hits   = np.frombuffer(overlap_hits.get_obj(),   dtype='u8').reshape(self.frame_shape)
        np_overlap_misses = np.frombuffer(overlap_misses.get_obj(), dtype='u8').reshape(self.frame_shape)
        
        nframes  = self.frames
        my_start = (nframes // self.nproc) * p
        my_end   = min((nframes // self.nproc) * (p+1), nframes)
        Nevents  = my_end - my_start 
        
        stime = time.time()
        f_vds = h5py.File(self.vds_file, 'r')
        for d in range(my_start, my_end):
            # read vds file (assume photon units)
            mask = f_vds[VDS_MASK_DATASET][d] == 0
             
            cids = f_vds['entry_1/cellId'][d]
            vals = mask * f_vds[VDS_DATASET][d]
            
            if self.is_hit[d] :
                np_powder_hits  += vals
                np_overlap_hits += mask
            
            if self.is_miss[d] :
                np_powder_misses  += vals
                np_overlap_misses += mask
            
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
    d = Path(out_fname_hits).resolve().parent
    assert(d.is_dir())
    
    # check that output directory exists
    d = Path(out_fname_misses).resolve().parent
    assert(d.is_dir())
    
    # load is_hit from events file
    with h5py.File(events_fname, 'r') as f:
        is_hit       = f['is_hit'][()]
        pulse_energy = f['pulse_energy'][()]
        hit_sigma    = f['hit_sigma'][()]
    
    # let's call a miss anything with a hitsigma bw +- 2
    sig_min = -2
    sig_max = 2
    is_miss = (hit_sigma > sig_min) * (hit_sigma < sig_max)
    
    m = pulse_energy > 0
    av_pulse_energy_hits   = np.mean(pulse_energy[m * is_hit])
    av_pulse_energy_misses = np.mean(pulse_energy[m * is_miss])
    
    calc_pow = Calc_powder(vds_file, is_hit, is_miss, nproc=args.nproc)
    
    powder_hits, powder_misses, overlap_hits, overlap_misses = calc_pow.run_frame()
    
    note_misses = f"""
    Misses are defined to be any frame with: {sig_min} < hit_sigma < {sig_max}
    """
    
    note_hits = f"""
    hits are defined by 'is_hit' in the events file
    """
    
    print(f'writing hits powder to {out_fname_hits}')
    with h5py.File(out_fname_hits, 'a') as f:
        utils.update_h5(f, 'powder_sum_hits', powder_hits, compression=True)
        utils.update_h5(f, 'overlap', overlap_hits, compression=True)
        utils.update_h5(f, 'is_hit', is_hit, compression=True)
        utils.update_h5(f, 'note', note_hits, compression=False)

    print(f'writing misses powder to {out_fname_misses}')
    with h5py.File(out_fname_misses, 'a') as f:
        utils.update_h5(f, 'powder_sum_misses', powder_hits, compression=True)
        utils.update_h5(f, 'overlap', overlap_misses, compression=True)
        utils.update_h5(f, 'is_miss', is_miss, compression=True)
        utils.update_h5(f, 'note', note_misses, compression=False)

