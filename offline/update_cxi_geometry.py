import argparse
from constants import PREFIX, DET_DIST, FRAME_SHAPE, VDS_DATASET, VDS_MASK_DATASET, SATURATION
import common

parser = argparse.ArgumentParser(description='write updated geometry to a cxi file')
parser.add_argument('run', type=int, nargs='+', help='Run number')

args = parser.parse_args()

for run in args.run :
    args.output_file   = PREFIX+'scratch/saved_hits/r%.4d_hits.cxi' %run
    args.geom_file     = common.get_geom(run)
    args.z             = DET_DIST

    import numpy as np
    import h5py
    import extra_geom



    # get pixle maps
    geom = extra_geom.AGIPD_1MGeometry.from_crystfel_geom(args.geom_file)
    xyz   = np.transpose(geom.get_pixel_positions(), (3, 0, 1, 2))
    x_pixel_size = geom.pixel_size
    y_pixel_size = geom.pixel_size
    pixel_area   = x_pixel_size * y_pixel_size
    xyz[2] = args.z

    print(f'updating {args.output_file} with {args.geom_file}' )

    with h5py.File(args.output_file, 'a') as f:
        key = '/entry_1/instrument_1/detector_1/xyz_map'
        f[key][:] = xyz

print('done')
