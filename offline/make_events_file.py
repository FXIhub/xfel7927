"""
add cellID, pulseID, trainID, detector distance to an events file, eg events/events_r0034.h5
"""

import os, sys
import argparse

import numpy as np
import h5py
import multiprocessing as mp
import glob
from pathlib import Path
from tqdm import tqdm
import utils
import shmemarray
import subprocess

PREFIX = os.environ["EXP_PREFIX"]
from constants import NMODULES, MODULE_SHAPE

def get_module_fnams(run, module):
    # get module file names
    run_dir = f'{PREFIX}/proc/r{run:>04}/'
    assert(Path(run_dir).is_dir())
    
    fnams = glob.glob(f'{run_dir}/CORR-*-AGIPD{module:>02}-S*.h5')
    
    print(f'found {len(fnams)} files for module {module}.')
    return fnams

def get_tid_cid_mapping_vds(vds):
    tc_to_vds_index = {}
    with h5py.File(vds) as f:
        tids = f['/entry_1/trainId'][()]
        cids = f['/entry_1/cellId'][:, 0]

    # initialise lookup 
    for t in np.unique(tids):
        tc_to_vds_index[t] = {}
    
    print(tids.shape, cids.shape)
    for i, (tid, cid) in enumerate(zip(tids, cids)):
        tc_to_vds_index[tid][cid] = i
    return tc_to_vds_index, tids.shape[0]


class Events():
    def __init__(self, module, fnams, mask, hit_mask, tc_to_vds_index, Nevents):
        self.fnams    = fnams
        self.module   = module
        self.mask     = mask
        self.hit_mask = hit_mask
        self.tc_to_vds_index = tc_to_vds_index
         
        self.data_src = f'/INSTRUMENT/SPB_DET_AGIPD1M-1/CORR/{module}CH0:output/image/data'

        N = 0
        # get total number of events
        for fnam in fnams:
            with h5py.File(fnam) as f:
                N += f[self.data_src].shape[0]
                self.frame_type = f[self.data_src].dtype
        
        self.events = Nevents
        
        print(f'found {N} events for module {module} in files and {self.events} events in VDS file for {len(self.fnams)} files')
        sys.stdout.flush()
        
        # check if we can skip hit finding mask
        if not np.any(hit_mask):
            self.skip_hit = True
        else :
            self.skip_hit = False
        
        self.lit      = shmemarray.empty((self.events,), dtype = np.uint64)
        self.lit_mask = shmemarray.empty((self.events,), dtype = np.uint64)
        self.lit[:]      = 0
        self.lit_mask[:] = 0
         
        self.counts      = shmemarray.empty((self.events,), dtype = np.uint64)
        self.counts_mask = shmemarray.empty((self.events,), dtype = np.uint64)
        self.counts[:]      = 0
        self.counts_mask[:] = 0
        
        self.cids = shmemarray.empty((self.events,), dtype = np.uint64)
        self.tids = shmemarray.empty((self.events,), dtype = np.uint64)
        self.cids[:] = 0
        self.tids[:] = 0
        

    def run(self):
        pool = mp.Pool(len(self.fnams))
        
        result_iter = pool.imap_unordered(
            self.run_worker, range(len(self.fnams)), chunksize = 1
        )
        
        for r in result_iter:
            pass
        
        self.finish()
        
        return self.lit, self.counts, self.lit_mask, self.counts_mask, self.events, self.cids, self.tids

    def finish(self):
        # I think this is how we go from shmemarray to np.array
        self.lit = np.array(self.lit)
        self.lit_mask = np.array(self.lit_mask)
        self.counts = np.array(self.counts)
        self.counts_mask = np.array(self.counts_mask)
        self.cids = np.array(self.cids)
        self.tids = np.array(self.tids)
    
    def run_worker(self, proc):
        fnam = self.fnams[proc]
        
        if proc == (len(self.fnams)-1) :
            disable = False
        else :
            disable = True
        
        data_src    = f'/INSTRUMENT/SPB_DET_AGIPD1M-1/CORR/{self.module}CH0:output/image/data'
        cellId_src  = f'/INSTRUMENT/SPB_DET_AGIPD1M-1/CORR/{self.module}CH0:output/image/cellId'
        trainId_src = f'/INSTRUMENT/SPB_DET_AGIPD1M-1/CORR/{self.module}CH0:output/image/trainId'
        processed   = 0
        skipped_events = 0
        with h5py.File(fnam) as f:
            data = f[data_src]
            cid  = f[cellId_src][()]
            tid  = f[trainId_src][()]
            
            frame = np.empty(MODULE_SHAPE, dtype = self.frame_type)
            
            for d in tqdm(range(data.shape[0]), disable = disable):
                try : 
                    i = self.tc_to_vds_index[tid[d]][cid[d]]
                except KeyError :
                    #print(f'skipping train {tid[d]} cell {cid[d]} (not present in vds file)')
                    skipped_events += 1
                    continue
                
                data.read_direct(frame, np.s_[d], np.s_[:])

                self.cids[i] = cid[d]
                self.tids[i] = tid[d]
                 
                # apply per-cell mask
                frame *= self.mask[cid[d]]
                
                self.lit[i]    = np.sum(frame>0)
                self.counts[i] = np.sum(frame)
                
                # apply hit_finding mask
                if self.skip_hit :
                    self.lit_mask[i]    = 0
                    self.counts_mask[i] = 0
                else :
                    frame *= self.hit_mask
                    self.lit_mask[i]    = np.sum(frame>0)
                    self.counts_mask[i] = np.sum(frame)
                
                processed += 1

        print('file', fnam, 'module', self.module, 'processed', processed, '/', self.events, 'events', 'skipped_events', skipped_events )
        sys.stdout.flush()
        return True


def main():
    PREFIX = os.environ["EXP_PREFIX"]
    parser = argparse.ArgumentParser(description='Lit pixel calculator')
    parser.add_argument('run', type=int, help='Run number')
    parser.add_argument('-o', '--out_folder', 
                        help='Path of output folder (default=%s/scratch/events/)'%PREFIX,
                        default=PREFIX+'/scratch/events/')
    parser.add_argument('--out_folder_powder', 
                        help='Path of output folder (default=%s/scratch/powder/)'%PREFIX,
                        default=PREFIX+'/scratch/powder/')
    parser.add_argument('-m', '--hit_finding_mask', 
                        help=f'pixel mask to apply before calculating litpixels and photon counts located in {PREFIX}scratch/det/. Mask must be 1 for good and 0 for bad, and the data must be in /entry_1/good_pixels)',
                        type=str, 
                        default='hit_finding_mask.h5' )
    args = parser.parse_args()
    
    # get autogenerated per-cell mask
    args.run_mask         = f'{PREFIX}scratch/det/r{args.run:>04}_mask.h5'
    args.hit_finding_mask = f'{PREFIX}scratch/det/{args.hit_finding_mask}'

    vds_file     = PREFIX+'scratch/vds/r%.4d.cxi' %args.run

    out_fname  = args.out_folder + os.path.splitext(os.path.basename(vds_file))[0] + '_events.h5'
    
    # check that run mask file exists
    assert(os.path.exists(args.run_mask))

    # check that vds file exists
    assert(os.path.exists(vds_file))
    
    # check that vds file points to corrected data
    result = subprocess.run(['h5ls', '-rv', vds_file], stdout=subprocess.PIPE)
    assert('CORR' in result.stdout.decode('utf-8'))

    print(f'loading hit finding mask {args.hit_finding_mask}')
    with h5py.File(args.hit_finding_mask) as f:
        hit_mask = f['entry_1/good_pixels'][()]
    
    # get lookup table for trainId, cellIds -> vds index
    tid_cid_mapping_vds, Nevents = get_tid_cid_mapping_vds(vds_file)
    
    litpixels     = []
    photon_counts = []
    cellIds       = []
    trainIds      = []
    modules       = []
    litpixels_mask     = []
    photon_counts_mask = []
    N_module = []

    
    for module in range(NMODULES):
        print(f'loading per-cell mask {args.run_mask} for module {module}')
        with h5py.File(args.run_mask) as f:
            cids = f['entry_1/cellIds'][()]
            mask = {}
            data = f['entry_1/good_pixels']
            for i, c in enumerate(cids) :
                mask[c] = data[i, module]
        
        fnams   = get_module_fnams(args.run, module)
        
        events = Events(module, fnams, mask, hit_mask[module], tid_cid_mapping_vds, Nevents)
        
        litpix, counts, litpix_mask, counts_mask, N, cids, tids = events.run()
        
        litpixels.append(litpix.copy())
        litpixels_mask.append(litpix_mask.copy())
        
        photon_counts.append(counts.copy())
        photon_counts_mask.append(counts_mask.copy())
        
        N_module.append(N)
        cellIds.append(cids.copy())
        trainIds.append(tids.copy())
        modules.append(module)

    # this is sure to fail, I just want to see when

    # check that we have the same number of events for each module
    assert(np.allclose(N_module[0], N_module))
    
    # check that cellIds agree between modules
    assert(np.allclose(cellIds[0], cellIds))

    # check that train agree between modules
    assert(np.allclose(trainIds[0], trainIds))
    
    # sum counts and litpix over modules
    litpixels_global     = np.sum(litpixels, axis=0)
    photon_counts_global = np.sum(photon_counts, axis=0)

    litpixels_mask_global     = np.sum(litpixels_mask, axis=0)
    photon_counts_mask_global = np.sum(photon_counts_mask, axis=0)
    
    print()
    print(f'writing output to {out_fname}')
    with h5py.File(out_fname, 'w') as f:
        utils.update_h5(f, 'total_intens', photon_counts_global, compression = True, chunks = photon_counts_global.shape)
        utils.update_h5(f, 'litpixels',    litpixels_global, compression = True, chunks = litpixels_global.shape)
        utils.update_h5(f, 'total_intens_mask', photon_counts_mask_global, compression = True, chunks = photon_counts_mask_global.shape)
        utils.update_h5(f, 'litpixels_mask',    litpixels_mask_global, compression = True, chunks = litpixels_mask_global.shape)
        utils.update_h5(f, 'cellId', cellIds[0], compression = True)
        utils.update_h5(f, 'trainId', trainIds[0], compression = True)
        utils.update_h5(f, 'modules', np.array(modules), compression = True)

if __name__ == '__main__':
    main()

