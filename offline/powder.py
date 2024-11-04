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
import sys

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
    
    return cellIds

def get_all_cell_ids(run):
    cs = []
    for module in range(NMODULES):
        fnams = get_module_fnams(run, module)
        cs.append(get_cell_ids(fnams, module))
    
    cellIds = np.unique(np.concatenate(cs))
    
    print()
    print(f'found {len(cellIds)} unique cellIds for all modules:')
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

def get_selection(select_events, run, out_folder):
    if select_events :
        assert(len(select_events) == 2)
        
        events   = f'{PREFIX}scratch/events/r{run:>04}_events.h5'
        
        # check that events file exists
        assert(os.path.exists(events))
        
        print(f'finding events in {events} where {select_events[0]} = {select_events[1]}')
        with h5py.File(events) as f:
            select = f[select_events[0]][()] == eval(select_events[1])
            tids = f['trainId'][()]
            cids = f['cellId'][()]
        
        print(f'found {np.sum(select)} events to process')
        
        # get lookup table for trainId, cellIds -> vds index
        tid_cid_mapping_vds = get_tid_cid_mapping_vds(tids, cids)

        select_tid_cid = {}
        for tid in tid_cid_mapping_vds.keys() :
            select_tid_cid[tid] = {}
            for cid in tid_cid_mapping_vds[tid].keys() :
                i = tid_cid_mapping_vds[tid][cid]
                if select[i] :
                    select_tid_cid[tid][cid] = True
                else :
                    select_tid_cid[tid][cid] = False
        
        out = f'{out_folder}/r{run:>04}_powder_{select_events[0]}_{select_events[1]}.h5'
    else :
        # this can't be pickled and so multiprocessing won't work
        # make a pretend nested dictionary
        #class A:
        #    def __getitem__(self, i):
        #        return True
        #class B:
        #    def __getitem__(self, i):
        #        return A()
        #select_tid_cid = B()
        select_tid_cid = None
 
        out = f'{out_folder}/r{run:>04}_powder.h5'

    return select_tid_cid, out
        


class Powdersum():
    def __init__(self, module, fnams, cellIds, select_tid_cid, nproc):
        # make inverse
        self.cid_to_index = {}
        for i, c in enumerate(cellIds) :
            self.cid_to_index[c] = i
        
        self.fnams  = fnams
        self.module = module
        
        self.powder_part = shmemarray.empty((len(fnams),) + (len(cellIds),) + MODULE_SHAPE, dtype = np.uint32)
        self.events_part = shmemarray.empty((len(fnams),) + (len(cellIds),)               , dtype = np.uint32)
        
        self.select_tid_cid = select_tid_cid
        self.nproc = nproc

    def run(self, fnams, module):
        self.fnams  = fnams
        self.module = module
        
        assert(len(fnams)<= self.powder_part.shape[0])

        if self.nproc is None :
            self.nproc = len(self.fnams)
        
        pool = mp.Pool(processes=self.nproc)
        
        result_iter = pool.imap_unordered(
            self.run_worker, range(len(self.fnams)), chunksize=1
        )
        
        for r in result_iter:
            pass

        self.finish()

        return self.powder, self.events
    
    def finish(self):
        F = len(self.fnams)
        self.powder = np.sum(self.powder_part[:F], axis=0)
        self.events = np.sum(self.events_part[:F], axis=0)
    
    def run_worker(self, proc):
        print(f'process {proc} starting')
        sys.stdout.flush() 
    
        fnam = self.fnams[proc]
        
        if proc == (len(self.fnams)-1) :
            disable = False
        else :
            disable = True
        
        self.powder_part[proc] = 0
        self.events_part[proc] = 0
         
        data_src    = f'/INSTRUMENT/SPB_DET_AGIPD1M-1/CORR/{self.module}CH0:output/image/data'
        cellId_src  = f'/INSTRUMENT/SPB_DET_AGIPD1M-1/CORR/{self.module}CH0:output/image/cellId'
        trainId_src = f'/INSTRUMENT/SPB_DET_AGIPD1M-1/CORR/{self.module}CH0:output/image/trainId'

        if self.select_tid_cid is not None :
            process_event = lambda t, c: self.select_tid_cid[t][c]
        else :
            process_event = lambda t, c: True
        
        with h5py.File(fnam) as f:
            data = f[data_src]
            cid  = f[cellId_src][()]
            tid  = f[trainId_src][()]
            
            frame = np.empty(MODULE_SHAPE, dtype = self.powder_part.dtype)
            
            for d in tqdm(range(data.shape[0]), disable = disable):
                t = tid[d]
                c = cid[d]
                try :
                    p = process_event(t, c) 
                except KeyError :
                    p = False

                if p :
                    data.read_direct(frame, np.s_[d], np.s_[:])
                    i = self.cid_to_index[c]
                    self.powder_part[proc, i] += frame
                    self.events_part[proc, i] += 1
        
        print(f'process {proc} finished')
        sys.stdout.flush() 
        return True


def main():
    parser = argparse.ArgumentParser(description='calculate powder per cell')
    parser.add_argument('run', type=int, help='Run number')
    parser.add_argument('-s', '--select_events', nargs='+', type = str,
                        help='The first argument is a dataset in the events file. The second argument is the flag to select e.g. "-s is_hit True" will calculate powder only for is_hit')
    parser.add_argument('-n', '--nproc', type = int,
                        help='number of processors to use (if None then will use 1 process for every module file')
    parser.add_argument('-o', '--out_folder', 
                        help='Path of output folder (default=%s/scratch/powder/)'%PREFIX,
                        default=PREFIX+'/scratch/powder/')
    args = parser.parse_args()
    #args.out = f'{args.out_folder}/r{args.run:>04}_powder.h5'
    
    powders  = []
    eventss  = []
    cellIdss = []
    modules  = []
    
    select_tid_cid, out = get_selection(args.select_events, args.run, args.out_folder)
    args.out = out
    print(f'output filename {args.out}')
    
    # get all cellIds
    cellIds = get_all_cell_ids(args.run)
    
    for module in range(NMODULES):
        fnams   = get_module_fnams(args.run, module)
        
        powdersum = Powdersum(module, fnams, cellIds, select_tid_cid, args.nproc)
        
        powder, events = powdersum.run(fnams, module)
        
        powders.append(powder.copy())
        eventss.append(events.copy())
        cellIdss.append(cellIds.copy())
        modules.append(module)
    
    powderg = np.array(powders)
    
    print()
    print(f'writing output to {args.out}')
    with h5py.File(args.out, 'w') as f:
        utils.update_h5(f, 'data', powderg, compression = True, chunks = powder.shape)
        utils.update_h5(f, 'cellIds', cellIds, compression = True)
        utils.update_h5(f, 'events', np.array(eventss), compression = True)
        utils.update_h5(f, 'modules', np.array(modules), compression = True)


if __name__ == '__main__':
    main()

