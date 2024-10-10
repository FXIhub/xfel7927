import argparse
from constants import PREFIX, DET_DIST, FRAME_SHAPE, VDS_DATASET, VDS_MASK_DATASET
import common

parser = argparse.ArgumentParser(description='size hits and write to cxi file and events file')
parser.add_argument('run', type=int, help='Run number')
parser.add_argument('--rmin', type=float, default=0., help='minimum radius from beam centre')
parser.add_argument('--rmax', type=float, default=0.02, help='maximum radius from beam centre')
parser.add_argument('--smin', type=float, default=1e-9, help='minimum particle diameter')
parser.add_argument('--smax', type=float, default=200e-9, help='maximum particle diameter')
parser.add_argument('--N_sizes', type=int, default=64, help='number of sizes to fit')
parser.add_argument('--N_angles', type=int, default=64, help='number of angles to fit')
parser.add_argument('--bins', type=int, default=4, help='number of pixels to bin along each dimension')

args = parser.parse_args()
    
args.output_file   = PREFIX+'scratch/saved_hits/r%.4d_hits.cxi' %args.run
args.events_file   = PREFIX+'scratch/events/r%.4d_events.h5' %args.run

import numpy as np
import h5py
#import extra_data
from tqdm import tqdm
import sys
import time
import utils
#import multiprocessing as mp
import scipy.constants as sc

from sizing_spheroid import Sizing

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


print(f'Sizing hits from {args.output_file}')

with h5py.File(args.output_file) as f:
    indices = list(range(f['entry_1/data_1/data'].shape[0]))
    mask    = f['/entry_1/instrument_1/detector_1/good_pixels'][()]
    xyz     = f['/entry_1/instrument_1/detector_1/xyz_map'][()]
    pixel_size = f['/entry_1/instrument_1/detector_1/x_pixel_size'][()]
    pixel_size = f['/entry_1/instrument_1/detector_1/x_pixel_size'][()]
    photon_energy = f['/entry_1/instrument_1/source_1/photon_energy'][()]
    
photon_energy = np.mean(photon_energy)

Nevents = len(indices)

#size = 64

# split frames over ranks
events_rank = np.linspace(0, Nevents, size+1).astype(int)

# generate lookup table
sizing = Sizing(
            mask, xyz, pixel_size, photon_energy, 
            args.rmin, args.rmax, args.N_angles,
            args.smin, args.smax, args.N_sizes, 
            polarisation='x', bin_size=args.bins)


def worker():
    my_indices = indices[events_rank[rank]: events_rank[rank+1]] 
    
    print(f'rank {rank} is processing indices {events_rank[rank]} to {events_rank[rank+1]}')
    sys.stdout.flush()
    
    if rank == 0 :
        it = tqdm(range(len(my_indices)), desc = f'Processing data from {args.output_file}')
    else :
        it = range(len(my_indices))
    
    out = []
    with h5py.File(args.output_file) as g:
        data = g[VDS_DATASET]
         
        # need to implement
        vds_index = g['/entry_1/vds_index'][()]
        
        for i in it:
            index = my_indices[i]
            
            av_diameter, long_axis_radius, short_axis_radius, theta_x, theta_z = sizing.size(data[index])
            
            out.append( [index, vds_index[index], short_axis_radius, long_axis_radius, theta_x, theta_z] )
            
    #b.put(out)
    return out

out = worker()

all_out = comm.gather(out, root=0)

if rank == 0 :
    #out = np.array(all_out).reshape(-1, 4)
    out = []
    for i in all_out :
        out += i
    
    # get short and long axis in order
    indices           = np.array([o[0] for o in out], dtype=int)
    i = np.argsort(indices)
    
    vds_indices       = np.array([o[1] for o in out], dtype=int)
    short_axis_radius = np.array([o[2] for o in out])
    long_axis_radius  = np.array([o[3] for o in out])
    theta_x           = np.array([o[4] for o in out])
    theta_z           = np.array([o[5] for o in out])
     
    short_axis_diameter = 2 * np.array(short_axis_radius)[i]
    long_axis_diameter  = 2 * np.array(long_axis_radius)[i]

    # output to cxi file 
    with h5py.File(args.output_file, 'r+') as f:
        assert(np.allclose(vds_indices[i], f['entry_1/vds_index'][()]))
        
        key = 'entry_1/sizing/short_axis_diameter'
        utils.update_h5(f, key, short_axis_diameter.astype(np.float32), compression = True)
        
        key = 'entry_1/sizing/long_axis_diameter'
        utils.update_h5(f, key, long_axis_diameter.astype(np.float32), compression = True)
        
        key = 'entry_1/sizing/theta_x'
        utils.update_h5(f, key, theta_x[i].astype(np.float32), compression = True)
        
        key = 'entry_1/sizing/theta_z'
        utils.update_h5(f, key, theta_z[i].astype(np.float32), compression = True)

    # output to events file 
    # get indicies
    with h5py.File(args.events_file, 'r+') as f:
        indices = np.arange(f['trainId'].shape[0])
        
        short_axis_diameter = np.zeros(indices.shape[0], dtype = np.float32)
        long_axis_diameter  = np.zeros(indices.shape[0], dtype = np.float32)
        
        short_axis_diameter[vds_indices] = 2 * short_axis_radius
        long_axis_diameter[vds_indices]  = 2 * long_axis_radius
        
        key = 'sizing/short_axis_diameter'
        utils.update_h5(f, key, short_axis_diameter, compression = True)
        
        key = 'sizing/long_axis_diameter'
        utils.update_h5(f, key, long_axis_diameter, compression = True)

    print('Done')
    sys.stdout.flush()

