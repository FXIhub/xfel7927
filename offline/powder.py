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



class Powdersum():
    def __init__(self, module, fnams, cellIds):
        # make inverse
        self.cid_to_index = {}
        for i, c in enumerate(cellIds) :
            self.cid_to_index[c] = i

        self.fnams = fnams
        self.module = module
        
        self.powder_part = shmemarray.empty((len(fnams),) + (len(cellIds),) + MODULE_SHAPE, dtype = np.uint32)
        self.events_part = shmemarray.empty((len(fnams),) + (len(cellIds),)               , dtype = np.uint32)

    def run(self):
        pool = mp.Pool(len(self.fnams))
        
        result_iter = pool.imap_unordered(
            self.run_worker, range(len(self.fnams))
        )
        
        for r in result_iter:
            pass

        self.finish()

        return self.powder, self.events

    def finish(self):
        self.powder = np.sum(self.powder_part, axis=0)
        self.events = np.sum(self.events_part, axis=0)
    
    def run_worker(self, proc):
        fnam = self.fnams[proc]

        if proc == (len(self.fnams)-1) :
            disable = False
        else :
            disable = True
        
        self.powder_part[proc] = 0
        self.events_part[proc] = 0
         
        data_src    = f'/INSTRUMENT/SPB_DET_AGIPD1M-1/CORR/{self.module}CH0:output/image/data'
        cellId_src  = f'/INSTRUMENT/SPB_DET_AGIPD1M-1/CORR/{self.module}CH0:output/image/cellId'
        with h5py.File(fnam) as f:
            data = f[data_src]
            cid  = f[cellId_src][()]
            
            frame = np.empty(MODULE_SHAPE, dtype = self.powder_part.dtype)
            
            for d in tqdm(range(data.shape[0]), disable = disable):
                data.read_direct(frame, np.s_[d], np.s_[:])
                i = self.cid_to_index[cid[d]]
                self.powder_part[proc, i] += frame
                self.events_part[proc, i] += 1
        return True


def main():
    parser = argparse.ArgumentParser(description='calculate powder per cell')
    parser.add_argument('run', type=int, help='Run number')
    parser.add_argument('-o', '--out_folder', 
                        help='Path of output folder (default=%s/scratch/powder/)'%PREFIX,
                        default=PREFIX+'/scratch/powder/')
    args = parser.parse_args()
    args.out = f'{args.out_folder}/r{args.run:>04}_powder.h5'
    
    powders  = []
    eventss  = []
    cellIdss = []
    modules  = []
    
    for module in range(NMODULES):
        fnams   = get_module_fnams(args.run, module)
        cellIds = get_cell_ids(fnams, module)
        
        powdersum = Powdersum(module, fnams, cellIds)
        
        powder, events = powdersum.run()
        
        powders.append(powder.copy())
        eventss.append(events.copy())
        cellIdss.append(cellIds.copy())
        modules.append(module)
    
    # make global powder allowing for missing cells
    cellIds = np.unique(np.concatenate(cellIdss))
    powderg = np.zeros((NMODULES,) + powder.shape, powder.dtype)
    for ci, c in enumerate(cellIds):
        for module in range(NMODULES):
            cj = np.where(cellIdss[module] == c)[0]
            if len(cj) == 1 :
                powderg[module, ci] = powders[module][cj[0]]
    
    print()
    print(f'writing output to {args.out}')
    with h5py.File(args.out, 'w') as f:
        utils.update_h5(f, 'data', powderg, compression = True, chunks = powder.shape)
        utils.update_h5(f, 'cellIds', cellIds, compression = True)
        utils.update_h5(f, 'events', np.array(eventss), compression = True)
        utils.update_h5(f, 'modules', np.array(modules), compression = True)


if __name__ == '__main__':
    main()

