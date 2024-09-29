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
from tqdm import tqdm

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

print(f'{np.sum(is_hit)} hits in extra data')

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

vds_events = trainId_vds.shape[0]

# some trains might be excluded in vds file, have to align
# assume trains are sorted
# and that all pulses are present (should check these)
tr = np.unique(trainId_vds)

C = len(np.unique(cellId_vds))

hit_out = np.zeros((vds_events,), dtype = bool)
mis_out = np.zeros((vds_events,), dtype = bool)
lit_out = np.zeros((vds_events,), dtype = litpixels.dtype)
pho_out = np.zeros((vds_events,), dtype = photon_counts.dtype)
tra_out = trainId_vds.copy()
cel_out = cellId_vds.copy()
pul_out = pulseId_vds.copy()

utils.put_a_in_b_by_train(is_hit, trainId, cellId, hit_out, trainId_vds, cellId_vds)
utils.put_a_in_b_by_train(is_miss, trainId, cellId, mis_out, trainId_vds, cellId_vds)
utils.put_a_in_b_by_train(litpixels, trainId, cellId, lit_out, trainId_vds, cellId_vds)
utils.put_a_in_b_by_train(photon_counts, trainId, cellId, pho_out, trainId_vds, cellId_vds)
utils.put_a_in_b_by_train(trainId, trainId, cellId, tra_out, trainId_vds, cellId_vds)
utils.put_a_in_b_by_train(cellId, trainId, cellId, cel_out, trainId_vds, cellId_vds)
utils.put_a_in_b_by_train(pulseId, trainId, cellId, pul_out, trainId_vds, cellId_vds)

assert(np.allclose(trainId_vds, tra_out))
assert(np.allclose(cellId_vds, cel_out))
assert(np.allclose(pulseId_vds, pul_out))

# write to events file
with h5py.File(out_fname, 'a') as f:
    utils.update_h5(f, 'is_hit', hit_out, compression=True)
    utils.update_h5(f, 'is_miss', mis_out, compression=True)
    
    utils.update_h5(f, 'total_intens', pho_out, compression=True)
    utils.update_h5(f, 'litpixels', lit_out, compression=True)
    
    utils.update_h5(f, 'trainId', trainId_vds, compression=True)
    utils.update_h5(f, 'cellId', cellId_vds, compression=True)
    utils.update_h5(f, 'pulseId', pulseId_vds, compression=True)
