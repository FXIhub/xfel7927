import numpy as np
import h5py
import pyqtgraph as pg
import os
import extra_geom
import common
import sizing_spheroid

PREFIX = os.environ['EXP_PREFIX']
run = 400

cxi_file = f'{PREFIX}/scratch/saved_hits/r{run:04}_hits.cxi'

with h5py.File(cxi_file) as f:
    xyz  = f['/entry_1/instrument_1/detector_1/xyz_map'][()]
    mask = f['/entry_1/instrument_1/detector_1/good_pixels'][()]
    wav  = np.mean(f['/entry_1/instrument_1/source_1/photon_wavelength'][()])
    pixel_size  = f['/entry_1/instrument_1/detector_1/x_pixel_size'][()]


# get q
qmap, solid, pol = sizing_spheroid.make_q_map_solid_angle_polarisation(xyz, pixel_size, wav)

q = np.sum(qmap**2, axis=0)**0.5

res = 1/q    

geom_file = common.get_geom(run)
geom = extra_geom.AGIPD_1MGeometry.from_crystfel_geom(geom_file)

# plot in nm
pg.show(1e9 * geom.position_modules(res)[0])


if __name__ == '__main__':
    # allow Control-C
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL) 
    pg.exec()


