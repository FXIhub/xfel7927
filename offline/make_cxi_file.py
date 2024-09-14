import argparse
from constants import PREFIX, DET_DIST, FRAME_SHAPE, VDS_DATASET, VDS_MASK_DATASET, SATURATION
import common

parser = argparse.ArgumentParser(description='Save hits in photon units to a cxi file')
parser.add_argument('run', type=int, help='Run number')
parser.add_argument('-s', '--sample_name',
                    help='name of sample',
                    type=str, default='DNA Pointer')

args = parser.parse_args()
    
args.output_file      = PREFIX+'scratch/saved_hits/r%.4d_hits.cxi' %args.run
args.vds_file      = PREFIX+'scratch/vds/r%.4d.cxi' %args.run
args.events_file   = PREFIX+'scratch/events/r%.4d_events.h5'%args.run
args.geom_file     = common.get_geom(args.run)
args.z             = DET_DIST

import numpy as np
import h5py
#import extra_data
from tqdm import tqdm
import sys
import time
import extra_geom
import scipy.constants as sc

import multiprocessing as mp

data_dtype = np.unit8

# get pixle maps
geom = extra_geom.AGIPD_1MGeometry.from_crystfel_geom(args.geom_file)
xyz   = np.transpose(geom.get_pixel_positions(), (3, 0, 1, 2))
x_pixel_size = 200e-6
y_pixel_size = 200e-6
pixel_area   = x_pixel_size * y_pixel_size
xyz[2] = args.z


"""
h5ls r0035_events.h5
/cellId                  Dataset {1055808, 16}
/is_hit                  Dataset {1055808}
/is_miss                 Dataset {1055808}
/litpixels               Dataset {1055808, 16}
/pulseId                 Dataset {1055808}
/pulse_energy            Dataset {1055808}
/total_intens            Dataset {1055808, 16}
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
    
    cellId_lit    = f['cellId'][:, 0]
    trainId_lit   = f['trainId'][()]
    pulseId_lit   = f['pulseId'][()]

    photons      = np.sum(f['total_intens'][()], axis=1)
    litpixels    = np.sum(f['litpixels'][()], axis=1)
    pulse_energy = f['pulse_energy'][()]
    wavelength   = f['wavelength'][()]

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
    entry_1['pulseId']  = pulseId_lit[indices]

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

    detector_1.create_dataset("mask", 
            shape=(Nevents,) + FRAME_SHAPE, 
            dtype=bool,
            chunks=(1,) + FRAME_SHAPE,
            compression='gzip',
            compression_opts=1,
            shuffle=True)
    
    # link /entry_1/data_1/data
    f["entry_1/data_1/data"] = h5py.SoftLink('/entry_1/instrument_1/detector_1/data')


#size = mp.cpu_count()
size = 16

# split frames over ranks
events_rank = np.linspace(0, Nevents, size+1).astype(int)

# load event ids for cross-checking
with h5py.File(args.vds_file) as g:
    cellId_vds = g['entry_1/cellId'][:, 0]
    trainId_vds = g['entry_1/trainId'][()]
    pulseId_vds = g['entry_1/pulseId'][()]

def worker(rank, lock):
    my_indices = indices[events_rank[rank]: events_rank[rank+1]] 
    
    print(f'rank {rank} is processing indices {events_rank[rank]} to {events_rank[rank+1]}')
    sys.stdout.flush()

    if rank == 0 :
        it = tqdm(range(len(my_indices)), desc = f'Processing data from {args.vds_file}')
    else :
        it = range(len(my_indices))

    frame_buf = np.empty((len(my_indices),) + FRAME_SHAPE, dtype=data_dtype)
    mask_buf  = np.empty((len(my_indices),) + FRAME_SHAPE, dtype=bool)

    with h5py.File(args.vds_file) as g:
        data = g[VDS_DATASET]
        mask = g[VDS_MASK_DATASET]
        
        for i in it:
            index = my_indices[i]
            
            frame_buf[i] = np.squeeze(data[index])
            sat = frame_buf[i] => SATURATION
            mask_buf[i]  = np.squeeze(mask[index] == 0)
            mask_buf[i, sat] = False
            
            # make sure we have the right event
            assert(cellId_vds[index] == cellId_lit[index])
            assert(trainId_vds[index] == trainId_lit[index])
            assert(pulseId_vds[index] == pulseId_lit[index])
            
            
    # take turns writing frame_buf to file 
    it = tqdm(range(len(my_indices)), desc = f'rank {rank} writing data to {args.output_file}')
    
    if lock.acquire() :
        with h5py.File(args.output_file, 'a') as f:
            for i in it :
                f['entry_1/instrument_1/detector_1/data'][events_rank[rank] + i] = frame_buf[i]
                f['entry_1/instrument_1/detector_1/mask'][events_rank[rank] + i] = mask_buf[i]
        
        sys.stdout.flush()
        lock.release()

lock = mp.Lock()
jobs = [mp.Process(target=worker, args=(m, lock)) for m in range(size)]
[j.start() for j in jobs]
[j.join() for j in jobs]
print('Done')
sys.stdout.flush()

