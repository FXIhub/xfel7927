'''Constants for beamtime p7076 at SPB'''

import os

#PREFIX = '/gpfs/exfel/exp/SQS/202302/p004098/scratch/'
PREFIX = os.environ("EXP_PREFIX")

VDS_DATASET      = '/entry_1/instrument_1/detector_1/data'
VDS_MASK_DATASET = '/entry_1/instrument_1/detector_1/mask'

NMODULES = 16
BAD_CELLIDS = [810]
DET_NAME = 'SPB_DET_AGIPD1M-1'


"""
MODULE_SHAPE = (512, 128)
NPULSES = 585
NCELLS = 800
#NPULSES = 400
#NCELLS = 400
CHUNK_SIZE = 1
#BAD_CELLIDS = list(range(0,17,2)) + list(range(484,499)) + [810]
ADU_PER_PHOTON = 5.

MASK_FNAME = PREFIX + '/geom/badpixel_background_mask_r0328.h5'
GAIN_FNAME = PREFIX + '/geom/gain_columns.npy'
CELLID_FNAME = PREFIX + '/geom/cellids_725.npy'

XGM_DA_NUM = 1
XGM_DATASET = "INSTRUMENT/SQS_DIAG1_XGMD/XGM/DOOCS:output/data/intensitySa3TD"
NBUNCHES_DA_NUM = 1
NBUNCHES_DATASET = '/CONTROL/SQS_DIAG1_XGMD/XGM/DOOCS/pulseEnergy/numberOfSa3BunchesActual/value'
DCONF_DA_NUM = 1
DCONF_DATASET = '/RUN/SQS_NQS_DSSC/FPGA/PPT_Q1/epcRegisterFilePath/value'
"""
