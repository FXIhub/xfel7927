# I need powder patterns
# got it use runs: 478 Co2 n he
#                  479 Co2 n
#                  480 he


import h5py
import numpy as np
import json
import os
import sizing_spheroid

PREFIX = os.environ['EXP_PREFIX']

runs = [478, 479, 480, 551]
labels = len(runs) * []

# make radial average
class Radial_average():
    def __init__(self, xyz, mask, radial_bin_size, min_rad, max_rad):
        self.r = (xyz[0]**2 + xyz[1]**2)**0.5
        self.radial_bin_size = radial_bin_size
        
        self.mask = mask.copy()
        self.mask[self.r < min_rad] = False 
        self.mask[self.r > max_rad] = False 
        
        # integer label for radial bin
        self.rl = np.rint(self.r[self.mask] / radial_bin_size).astype(int).ravel()
        
        # just to speed things up a little
        self.rl = np.ascontiguousarray(self.rl)
        
        # number of pixels contributing to each bin 
        self.rcounts = np.bincount(self.rl.ravel())
        
        # for reference
        self.rs = np.arange(np.max(self.rl)+1) * self.radial_bin_size
        
    def make_rad_av(self, ar):
        rsum = np.bincount(self.rl, ar[self.mask])
        
        # normalise, might not want to do this
        rsum /= np.clip(self.rcounts, 1e-20, None)
        return rsum

# load Sample name from run_table
fnam = f'{PREFIX}/scratch/log/run_table.json'

with open(fnam, 'r') as f:
    run_table = json.load(f)


labels = [run_table[str(run)]['Sample'] for run in runs]

# load powders 
fnams = [f'{PREFIX}/usr/Shared/amorgan/geomtools/nb/powder_sum_p007927_r{run:04}.h5' for run in runs]

powders = []
for fnam in fnams:
    with h5py.File(fnam) as f:
        powders.append(f['/powderSum/image/mean'][()])

# get radial profiles 
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
#res = 1/q    

# make radial average
r = Radial_average(xyz, mask, 4 * 200e-6, 0, 1000)

radial_profiles = [r.make_rad_av(p) for p in powders]
q_1d = r.make_rad_av(q)

# get pulse energies
av_pulse_energy = []
for run in runs:
    fnam_events = f'{PREFIX}/scratch/events/r{run:04}_events.h5'
    with h5py.File(fnam_events) as f:
        av_pulse_energy.append(np.mean(f['/pulse_energy'][()]))


import matplotlib.pyplot as plt
fig, ax = plt.subplots()
fig.set_size_inches(10, 3)
fig.set_tight_layout(True)


for run, label, rad, ave in zip(runs, labels, radial_profiles, av_pulse_energy) :
    print(q_1d.shape, rad.shape)
    ax.plot(1e-9 *q_1d, rad/(1e3 * ave), label= f'{label} run: {run}')

ax.legend(markerscale=5)
ax.set_yscale('log')
ax.spines[['right', 'top']].set_visible(False)
ax.set_ylabel('photon counts / mJ')
ax.set_xlabel('q (nm)^-1')
ax.grid(visible=True, which='both', alpha = 0.3)
plt.show()
