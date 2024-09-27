"""
add cellID, pulseID, trainID, detector distance to an events file, eg events/events_r0034.h5
"""

# doing this first speeds up running
import os
import argparse

from constants import EXP_ID, PREFIX

parser = argparse.ArgumentParser(description='Lit pixel calculator')
parser.add_argument('run', type=int, help='Run number')
parser.add_argument('-o', '--out_folder', 
                    help='Path of output folder (default=%s/scratch/events/)'%PREFIX,
                    default=PREFIX+'/scratch/events/')
args = parser.parse_args()
vds_file   = PREFIX+'scratch/vds/r%.4d.cxi' %args.run
out_fname  = args.out_folder + os.path.splitext(os.path.basename(vds_file))[0] + '_events.h5'

from extra_data import open_run
import h5py
import numpy as np
import utils

run = open_run(EXP_ID, args.run, data='all')

# hitfindig results
# run['SPB_DET_AGIPD1M-1/REDU/SPI_HITFINDER:output'].keys()
# data.hitFlag - boolean indicates hits
# data.missFlag - boolean indicates selected misses for background characterization
# data.dataFlag = data.hitFlag | data.missFlag
# data.hitscore - hitscore
# data.pulseId - pulse Id
# data.trainId - train Id
hit = run['SPB_DET_AGIPD1M-1/REDU/SPI_HITFINDER:output']
is_hit  = hit['data.hitFlag'].ndarray()
is_miss = hit['data.missFlag'].ndarray()
trainId = hit['data.trainId'].ndarray()
pulseId = hit['data.pulseId'].ndarray()

# hitfinder configuration
# run['SPB_DET_AGIPD1M-1/REDU/SPI_HITFINDER'].keys()

# litpixels
# run['SPB_DET_AGIPD1M-1/CORR/0CH0:output'].keys()
a = run['SPB_DET_AGIPD1M-1/CORR/0CH0:output']
cellId        = a['litpx.cellId'].ndarray()
litpixels     = a['litpx.litPixels'].ndarray()
photon_counts = a['litpx.totalIntensity'].ndarray()

# check that id's match vds
with h5py.File(vds_file, 'r') as f_vds:
    trainId_vds = f_vds['entry_1/trainId'][:]
    cellId_vds  = f_vds['entry_1/cellId'][:, 0]
    pulseId_vds = f_vds['entry_1/pulseId'][:]


assert(np.allclose(trainId_vds, trainId))
assert(np.allclose(cellId_vds, cellId))
assert(np.allclose(pulseId_vds, pulseId))

# write to events file
with h5py.File(out_fname, 'a') as f:
    utils.update_h5(f, 'is_hit', is_hit.astype(bool), compression=True)
    utils.update_h5(f, 'is_miss', is_miss.astype(bool), compression=True)
    
    utils.update_h5(f, 'total_intens', np.array(photon_counts), compression=True)
    utils.update_h5(f, 'litpixels', np.array(litpixels), compression=True)
    
    utils.update_h5(f, 'trainId', trainId, compression=True)
    utils.update_h5(f, 'cellId', cellId, compression=True)
    utils.update_h5(f, 'pulseId', pulseId, compression=True)
