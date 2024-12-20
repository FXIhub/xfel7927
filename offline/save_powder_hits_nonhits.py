# use multiprocessing to calculate powder pattern
# per module per cell 
# use Egor's shmem thing (looks nice)
# use mp.Pool

import os
import argparse

import numpy as np
import h5py
import multiprocessing as mp
import glob
from pathlib import Path
from tqdm import tqdm
import utils
import shmemarray

PREFIX = os.environ["EXP_PREFIX"]
from constants import NMODULES, MODULE_SHAPE

def get_module_fnams(run, module):
    # get module file names
    run_dir = f'{PREFIX}/proc/r{run:>04}/'
    assert(Path(run_dir).is_dir())
    
    fnams = glob.glob(f'{run_dir}/CORR-*-AGIPD{module:>02}-S*.h5')
    
    print(f'found {len(fnams)} files for module {module}.')
    return fnams

def get_cell_ids(fnams, module):
    src  = f'/INSTRUMENT/SPB_DET_AGIPD1M-1/CORR/{module}CH0:output/image/cellId'
    cs = []
    for fnam in fnams:
        with h5py.File(fnam) as f:
            cs.append(np.unique(f[src][()]))
    cellIds = np.unique(np.concatenate(cs))

    print()
    print(f'found {len(cellIds)} unique cellIds:')
    print(cellIds)
    return cellIds

def get_tid_cid_mapping_vds(tids, cids):
    tc_to_vds_index = {}

    # initialise lookup 
    for t in np.unique(tids):
        tc_to_vds_index[t] = {}
    
    for i, (tid, cid) in enumerate(zip(tids, cids)):
        tc_to_vds_index[tid][cid] = i
    return tc_to_vds_index



class Powdersum():
    def __init__(self, module, fnams, mask, is_hit, tid_cid_mapping_vds):
        self.fnams = fnams
        self.module = module
        
        self.powder_hit_part    = shmemarray.empty((len(fnams),) + MODULE_SHAPE, dtype = np.uint32)
        self.powder_nonhit_part = shmemarray.empty((len(fnams),) + MODULE_SHAPE, dtype = np.uint32)
        
        self.overlap_hit_part    = shmemarray.empty((len(fnams),) + MODULE_SHAPE, dtype = np.uint32)
        self.overlap_nonhit_part = shmemarray.empty((len(fnams),) + MODULE_SHAPE, dtype = np.uint32)

        self.Nhits_part    = shmemarray.empty((len(fnams),), dtype = np.uint32)
        self.Nnonhits_part = shmemarray.empty((len(fnams),), dtype = np.uint32)
        
        self.mask = mask
        self.is_hit = is_hit
        self.tid_cid_mapping_vds = tid_cid_mapping_vds

    def run(self):
        pool = mp.Pool(len(self.fnams))
        
        result_iter = pool.imap_unordered(
            self.run_worker, range(len(self.fnams))
        )
        
        for r in result_iter:
            pass
        
        self.finish()
        
        return self.powder_hit, self.powder_nonhit, self.overlap_hit, self.overlap_nonhit, self.Nhits, self.Nnonhits

    def finish(self):
        self.powder_hit     = np.sum(self.powder_hit_part, axis=0)
        self.powder_nonhit  = np.sum(self.powder_nonhit_part, axis=0)
        self.overlap_hit    = np.sum(self.overlap_hit_part, axis=0)
        self.overlap_nonhit = np.sum(self.overlap_nonhit_part, axis=0)
        self.Nhits          = np.sum(self.Nhits_part, axis=0)
        self.Nnonhits       = np.sum(self.Nnonhits_part, axis=0)

    def run_worker(self, proc):
        fnam = self.fnams[proc]

        if proc == (len(self.fnams)-1) :
            disable = False
        else :
            disable = True
        
        self.powder_hit_part[proc] = 0
        self.powder_nonhit_part[proc] = 0
        self.overlap_hit_part[proc] = 0
        self.overlap_nonhit_part[proc] = 0
        self.Nhits_part[proc] = 0
        self.Nnonhits_part[proc] = 0
         
        data_src    = f'/INSTRUMENT/SPB_DET_AGIPD1M-1/CORR/{self.module}CH0:output/image/data'
        cellId_src  = f'/INSTRUMENT/SPB_DET_AGIPD1M-1/CORR/{self.module}CH0:output/image/cellId'
        trainId_src = f'/INSTRUMENT/SPB_DET_AGIPD1M-1/CORR/{self.module}CH0:output/image/trainId'
        with h5py.File(fnam) as f:
            data = f[data_src]
            cid  = f[cellId_src][()]
            tid  = f[trainId_src][()]
            
            frame = np.empty(MODULE_SHAPE, dtype = self.powder_hit_part.dtype)
            
            for d in tqdm(range(data.shape[0]), disable = disable):
                data.read_direct(frame, np.s_[d], np.s_[:])
                frame *= self.mask[cid[d]]
                is_hit = self.is_hit[self.tid_cid_mapping_vds[tid[d]][cid[d]]]
                if is_hit :
                    self.powder_hit_part[proc]  += frame
                    self.overlap_hit_part[proc] += self.mask[cid[d]]
                    self.Nhits_part[proc]       += 1
                else :
                    self.powder_nonhit_part[proc]  += frame
                    self.overlap_nonhit_part[proc] += self.mask[cid[d]]
                    self.Nnonhits_part[proc]       += 1
        return True


def main():
    parser = argparse.ArgumentParser(description='calculate per-pixel powder patter for all hits and all non-hits')
    parser.add_argument('run', type=int, help='Run number')
    parser.add_argument('-o', '--out_folder', 
                        help='Path of output folder (default=%s/scratch/powder/)'%PREFIX,
                        default=PREFIX+'/scratch/powder/')
    args = parser.parse_args()
    args.out  = f'{args.out_folder}/r{args.run:>04}_powder_hits.h5'

    args.run_mask = f'{PREFIX}scratch/det/r{args.run:>04}_mask.h5'
    args.events   = f'{PREFIX}scratch/events/r{args.run:>04}_events.h5'

    # check that run mask and events file exists
    assert(os.path.exists(args.run_mask))
    assert(os.path.exists(args.events))

    with h5py.File(args.events) as f:
        is_hit = f['is_hit'][()]
        tids = f['trainId'][()]
        cids = f['cellId'][()]
    
    # get lookup table for trainId, cellIds -> vds index
    tid_cid_mapping_vds = get_tid_cid_mapping_vds(tids, cids)
    
    powders_hits     = []
    powders_nonhits  = []
    overlaps_hits     = []
    overlaps_nonhits     = []
    events_hits  = []
    events_nonhits  = []
    
    modules  = []

    for module in range(NMODULES):
        print(f'loading per-cell mask {args.run_mask} for module {module}')
        with h5py.File(args.run_mask) as f:
            cids = f['entry_1/cellIds'][()]
            mask = {}
            data = f['entry_1/good_pixels']
            for i, c in enumerate(cids) :
                mask[c] = data[i, module]

        fnams   = get_module_fnams(args.run, module)
        
        powdersum = Powdersum(module, fnams, mask, is_hit, tid_cid_mapping_vds)
        
        powder_hits, powder_nonhits, overlap_hits, overlap_non_hits, Nhits, Nnonhits = powdersum.run()
        
        powders_hits.append(powder_hits.copy())
        powders_nonhits.append(powder_nonhits.copy())

        overlaps_hits.append(overlap_hits.copy())
        overlaps_nonhits.append(overlap_non_hits.copy())
        
        events_hits.append(Nhits)
        events_nonhits.append(Nnonhits)
        modules.append(module)
    
    p = np.array(powders_hits)
    o = np.clip(np.array(overlaps_hits), 1, None)
    mean_hits = p / o
    
    p = np.array(powders_nonhits)
    o = np.clip(np.array(overlaps_nonhits), 1, None)
    mean_nonhits = p / o
    
    p = np.array(powders_hits) + np.array(powders_nonhits)
    o = np.clip(np.array(overlaps_hits) + np.array(overlaps_nonhits), 1, None)
    mean_total = p / o

    events_hits = np.array(events_hits)
    events_nonhits = np.array(events_nonhits)
    events = events_hits + events_nonhits

    note = f"""
    hits are defined by 'is_hit' in the events file {args.events}
    per-cell mask is defined by 'entry_1/good_pixels' in the file {args.run_mask}
    """
    
    print()
    print(f'writing output to {args.out}')
    with h5py.File(args.out, 'w') as f:
        utils.update_h5(f, 'powder_mean_hits',    mean_hits,    compression = True, chunks = mean_hits.shape)
        utils.update_h5(f, 'powder_mean_nonhits', mean_nonhits, compression = True, chunks = mean_hits.shape)
        utils.update_h5(f, 'powder_mean_total',   mean_total,   compression = True, chunks = mean_hits.shape)
        
        utils.update_h5(f, 'events', events, compression = True)
        utils.update_h5(f, 'events_hits', events_hits, compression = True)
        utils.update_h5(f, 'events_nonhits', events_nonhits, compression = True)
        utils.update_h5(f, 'modules', np.array(modules), compression = True)
         
        utils.update_h5(f, 'note', note, compression=False)

if __name__ == '__main__':
    main()

