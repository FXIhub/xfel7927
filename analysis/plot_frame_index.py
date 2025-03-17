import argparse
import extra_geom
import h5py

from constants import PREFIX, VDS_DATASET 


def parse_cmdline_args():
    parser = argparse.ArgumentParser(description='view frames from saved VDS files')
    parser.add_argument('run', type=int, help="run number of the VDS file.")
    parser.add_argument('index', type=int, help="frame index in file")
    
    parser.add_argument('--vmin', type=float, default = None, help="minimum value")
    parser.add_argument('--vmax', type=float, default = None, help="maximum value")
    
    args = parser.parse_args()
    args.vds_file    = PREFIX+'scratch/vds/r%.4d.cxi' % args.run
    return args

# for (much faster) local viewing
geom_fnam=f'../geom/r0600.geom'

args = parse_cmdline_args()

# get frame 
with h5py.File(args.vds_file) as f:
    frame = f[VDS_DATASET][args.index]

# display
geom = extra_geom.AGIPD_1MGeometry.from_crystfel_geom(geom_fnam)

import matplotlib.pyplot as plt
geom.plot_data(frame, vmin=args.vmin, vmax=args.vmax)
plt.show()
