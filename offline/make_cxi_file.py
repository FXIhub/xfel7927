import argparse
from constants import PREFIX, DET_DIST, FRAME_SHAPE, VDS_DATASET, VDS_MASK_DATASET, SATURATION
import common

parser = argparse.ArgumentParser(description='Save hits in photon units to a cxi file')
parser.add_argument('run', type=int, help='Run number')
#parser.add_argument('-s', '--sample_name',
#                    help='name of sample',
#                    type=str, default='DNA Pointer')
parser.add_argument('-m', '--mask', type=str, help=f'By default per-cell masks in {PREFIX}/scratch/det/r<runno>_mask.h5 are applied. Add a filename of global good pixels mask, located in {PREFIX}/scratch/det/')
parser.add_argument('-n', '--nproc', type=int, default=1, help=f'number of processes to use')

args = parser.parse_args()


# get sample name from run_table
# load run table
import os
import json
default_run_table = f'{PREFIX}scratch/log/run_table.json'
print(f'getting sample name from run table')
run_table = json.load(open(default_run_table, 'r'))
args.sample_name = None
for r, v in run_table.items():
    if not isinstance(v, dict):
        continue
    
    if v['Run number'] == args.run :
        args.sample_name = v['Sample']

assert(args.sample_name is not None)
print('sample name:', args.sample_name)


    
args.output_file   = PREFIX+'scratch/saved_hits/r%.4d_hits.cxi' %args.run
args.vds_file      = PREFIX+'scratch/vds/r%.4d.cxi' %args.run
args.events_file   = PREFIX+'scratch/events/r%.4d_events.h5'%args.run
if args.mask :
    args.mask_file     = f'{PREFIX}scratch/det/{args.mask}'
else:
    args.mask_file     = None
args.geom_file     = common.get_geom(args.run)
args.z             = DET_DIST

args.run_mask         = f'{PREFIX}scratch/det/r{args.run:>04}_mask.h5'
# check that run mask file exists
assert(os.path.exists(args.run_mask))

import numpy as np
import h5py
#import extra_data
from tqdm import tqdm
import sys
import time
import extra_geom
import scipy.constants as sc

import multiprocessing as mp

data_dtype = np.uint8

# get pixle maps
geom = extra_geom.AGIPD_1MGeometry.from_crystfel_geom(args.geom_file)
xyz   = np.transpose(geom.get_pixel_positions(), (3, 0, 1, 2))
x_pixel_size = 200e-6
y_pixel_size = 200e-6
pixel_area   = x_pixel_size * y_pixel_size
xyz[2] = args.z

# load per-cell mask
print(f'loading per-cell mask {args.run_mask}')
with h5py.File(args.run_mask) as f:
    cids = f['entry_1/cellIds'][()]
    mask = {}
    data = f['entry_1/good_pixels']
    for i, c in enumerate(cids) :
        mask[c] = data[i]

# add pixels for which all cells are masked to good_pixels
temp = np.ones(xyz[0].shape, dtype = int)
for c in mask.keys():
    temp += mask[c].astype(int)    

good_pixels = temp > 0

# get mask
if args.mask_file :
    with h5py.File(args.mask_file) as f:
        good_pixels *= f['entry_1/good_pixels'][()]

print('{round(100 * np.sum(~good_pixels) / good_pixels.size, 2)}% masked pixels in good_pixels')
sys.exit()
"""
h5ls r0035_events.h5
/cellId                  Dataset {1055808, 16}
/is_hit                  Dataset {1055808}
/is_miss                 Dataset {1055808}
/litpixels               Dataset {1055808}
/hitscore               Dataset {1055808}
/pulseId                 Dataset {1055808}
/pulse_energy            Dataset {1055808}
/total_intens            Dataset {1055808}
/trainId                 Dataset {1055808}
/wavelength              Dataset {1055808}
"""

print(f'loading hit selection from {args.events_file}')
with h5py.File(args.events_file) as f:
    # indices for definite hits
    m = f['is_hit'][()]
    indices = np.where(m)[0]
    
    # indices for definite miss
    m = f['is_miss'][()]
    indices_miss = np.where(m)[0]
        
    if len(f['cellId'].shape) == 2 :
        cellId_lit    = f['cellId'][:, 0]
    elif len(f['cellId'].shape) == 1 :
        cellId_lit    = f['cellId'][:]
    else :
        raise ValueError(f'unknown cellId shape in {args.events_file}')
    
    trainId_lit   = f['trainId'][()]
    # hack
    #pulseId_lit   = trainId_lit   #f['pulseId'][()]

    photons      = f['total_intens'][()]
    litpixels    = f['litpixels'][()]
    pulse_energy = f['pulse_energy'][()]
    wavelength   = f['wavelength'][()]

    if 'hit_score' in f :
        hitscore    = f['hit_score'][()]
    else :
        hitscore    = None

    if 'hit_sigma' in f :
        hit_sigma    = f['hit_sigma'][()]
    else :
        hit_sigma    = None

print(wavelength.shape)
photon_energy = sc.h * sc.c / wavelength
    
Nevents = len(indices)

print(f'found {Nevents} events labeled as hit')
sys.stdout.flush()

# entry_1/
#   name = "EUXFEL 2023 P3004 {run}"
#   experiment_identifier = index 
#   trainID
#   cellID
#   pulseID
#   start_time = {start time}
#   sample_1/
#       name "DNA Pointer"
#   instrument_1
#       source_1/
#           photon_energy
#           photon_wavelength
#           pulse_energy
#       detector_1/
#           data 
#           mask 
#           good_pixels 
#           xyz_map 
#           x_pixel_size 
#           y_pixel_size 
#           pixel_area 
#           background 
#           photon counts 
#           lit pixels 

print('Initialising cxi file ', args.output_file)
sys.stdout.flush()

with h5py.File(args.output_file, 'w') as f:
    entry_1 = f.create_group("entry_1")
    entry_1["name"] = f'P7927_DNA_Scaffold_Nanofocus_SPI_SPB_Sep_2024 run {args.run}'
    
    # copy sample name
    sample = entry_1.create_group("sample_1")
    sample['name'] = args.sample_name

    entry_1['trainId']  = trainId_lit[indices]
    entry_1['cellId']   = cellId_lit[indices]
    #entry_1['pulseId']  = pulseId_lit[indices]
    entry_1['vds_index']  = indices

    # copy instrument name
    instrument_1 = entry_1.create_group("instrument_1")
    instrument_1["name"] = 'SPB'
    
    entry_1['experiment_identifier'] = indices
    
    source_1 = instrument_1.create_group("source_1")
    source_1.create_dataset("photon_energy",     data = photon_energy[indices].astype(np.float32), compression='gzip', compression_opts=1, shuffle=True, chunks = True)
    source_1.create_dataset("photon_wavelength", data = wavelength[indices].astype(np.float32), compression='gzip', compression_opts=1, shuffle=True, chunks = True)
    source_1.create_dataset("pulse_energy",      data = pulse_energy[indices].astype(np.float32), compression='gzip', compression_opts=1, shuffle=True, chunks = True)
    source_1['photon_energy'].attrs['axes'] = "experiment_identifier"
    source_1['photon_wavelength'].attrs['axes'] = "experiment_identifier"
    source_1['pulse_energy'].attrs['axes'] = "experiment_identifier"
    
    data_1 = instrument_1.create_group("data_1")
    detector_1 = instrument_1.create_group("detector_1")

    detector_1['x_pixel_size'] = x_pixel_size
    detector_1['y_pixel_size'] = y_pixel_size
    detector_1['pixel_area'] = x_pixel_size * y_pixel_size
    
    detector_1.create_dataset("photon_counts",  data = photons[indices].astype(np.float32), compression='gzip', compression_opts=1, shuffle=True, chunks = True)
    detector_1.create_dataset("lit_pixels",   data = litpixels[indices].astype(np.float32), compression='gzip', compression_opts=1, shuffle=True, chunks = True)
    detector_1['photon_counts'].attrs['axes'] = "experiment_identifier"
    detector_1['lit_pixels'].attrs['axes'] = "experiment_identifier"

    if hitscore is not None :
        detector_1.create_dataset("hit_score",   data = hitscore[indices].astype(np.float32), compression='gzip', compression_opts=1, shuffle=True, chunks = True)
        detector_1['hit_score'].attrs['axes'] = "experiment_identifier"

    if hit_sigma is not None :
        detector_1.create_dataset("hit_sigma",   data = hit_sigma[indices].astype(np.float32), compression='gzip', compression_opts=1, shuffle=True, chunks = True)
        detector_1['hit_sigma'].attrs['axes'] = "experiment_identifier"
    
    # write pixel map
    detector_1.create_dataset('xyz_map', data = xyz, compression='gzip', compression_opts=1, shuffle = True, dtype = np.float32)
    
    # it's fine to compress if we are using pickle later
    detector_1.create_dataset("data", 
            shape=(Nevents,) + FRAME_SHAPE, 
            dtype=data_dtype,
            chunks=(1,) + FRAME_SHAPE,
            compression='gzip',
            compression_opts=1,
            shuffle=True)

    #detector_1.create_dataset("powder", 
    #        shape=FRAME_SHAPE, 
    #        dtype=np.uint64,
    #        chunks=FRAME_SHAPE,
    #        compression='gzip',
    #        compression_opts=1,
    #        shuffle=True)

    #detector_1.create_dataset("powder_overlap", 
    #        shape=FRAME_SHAPE, 
    #        dtype=np.uint64,
    #        chunks=FRAME_SHAPE,
    #        compression='gzip',
    #        compression_opts=1,
    #        shuffle=True,
    #        fillvalue = 0)

    detector_1.create_dataset("good_pixels", 
            data=good_pixels, 
            chunks=FRAME_SHAPE,
            compression='gzip',
            compression_opts=1,
            shuffle=True,
            fillvalue = 0)
    print('writing mask with {round(100 * np.sum(~good_pixels) / good_pixels.size, 2)}% of pixels masked')
    
    # link /entry_1/data_1/data
    f["entry_1/data_1/data"] = h5py.SoftLink('/entry_1/instrument_1/detector_1/data')

#size = mp.cpu_count()
size = args.nproc

# split frames over ranks
events_rank = np.linspace(0, Nevents, size+1).astype(int)

# load event ids for cross-checking
with h5py.File(args.vds_file) as g:
    cellId_vds = g['entry_1/cellId'][:, 0]
    trainId_vds = g['entry_1/trainId'][()]
    #pulseId_vds = g['entry_1/pulseId'][()]

def worker(rank, lock):
    my_indices = indices[events_rank[rank]: events_rank[rank+1]] 
    
    print(f'rank {rank} is processing indices {events_rank[rank]} to {events_rank[rank+1]}')
    sys.stdout.flush()
    
    if rank == 0 :
        it = tqdm(range(len(my_indices)), desc = f'Processing data from {args.vds_file}')
    else :
        it = range(len(my_indices))
    #powder    = np.zeros(FRAME_SHAPE, dtype = np.uint64)
    #overlap   = np.zeros(FRAME_SHAPE, dtype = np.uint64)
    
    frame_buf = np.empty((len(my_indices),) + FRAME_SHAPE, dtype=data_dtype)

    with h5py.File(args.vds_file) as g:
        data = g[VDS_DATASET]
        #mask = g[VDS_MASK_DATASET]
        
        for i in it:
            index = my_indices[i]
            cid   = cellId_vds[index]
            
            frame_buf[i] = np.squeeze(data[index])
            #m            = np.squeeze(mask[index] == 0)
              
            # apply per-frame mask
            frame_buf[i] *= mask[cid]  
            
            # add to powder
            #powder  += frame_buf[i]
            #overlap += m
            
            # make sure we have the right event
            assert(cellId_vds[index] == cellId_lit[index])
            assert(trainId_vds[index] == trainId_lit[index])
            #assert(pulseId_vds[index] == pulseId_lit[index])
            
    # take turns writing frame_buf to file 
    it = tqdm(range(len(my_indices)), desc = f'rank {rank} writing data to {args.output_file}')
    
    if lock.acquire() :
        with h5py.File(args.output_file, 'a') as f:
            for i in it :
                f['entry_1/instrument_1/detector_1/data'][events_rank[rank] + i] = np.clip(frame_buf[i], 0, np.iinfo(frame_buf.dtype).max)
            
            # update powder
            #powder_file  = f['entry_1/instrument_1/detector_1/powder'][()]
            #powder_file += powder
            #f['entry_1/instrument_1/detector_1/powder'][:] = powder_file
                
            #overlap_file  = f['entry_1/instrument_1/detector_1/powder_overlap'][()]
            #overlap_file += overlap
            #f['entry_1/instrument_1/detector_1/powder_overlap'][:] = overlap_file
        
        sys.stdout.flush()
        lock.release()

lock = mp.Lock()
jobs = [mp.Process(target=worker, args=(m, lock)) for m in range(size)]
[j.start() for j in jobs]
[j.join() for j in jobs]
print('Done')
sys.stdout.flush()

