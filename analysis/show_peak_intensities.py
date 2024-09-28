import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
import glob

PREFIX = os.environ["EXP_PREFIX"]

sample = 'cube'
size_min = 20e-9
size_max = 50e-9

# load hit_finding mask
fnam = f'{PREFIX}/scratch/det/hit_finding_mask.h5'
with h5py.File(fnam) as f:
    hit_mask = f['entry_1/good_pixels'][()]

# loop over runs in saved hits
fnams = sorted(glob.glob(f'{PREFIX}/scratch/saved_hits/*.cxi'))


runs = []
photons = []
for fnam in fnams:
    run = int(fnam[-13:-9])
    with h5py.File(fnam) as f:
        if '/entry_1/sizing' in f:
            sizes  = f['/entry_1/sizing/long_axis_diameter'][()] 
            sizes += f['/entry_1/sizing/short_axis_diameter'][()] 
            sizes /= 2
               
            mask  = f['/entry_1/instrument_1/detector_1/good_pixels'][()] 
            
            data = f['/entry_1/instrument_1/detector_1/data']
            
            ds = np.where( (sizes > size_min) * (sizes < size_max) )[0]
            
            if len(ds) > 0 :
                for d in ds :
                    photons.append(np.sum(data[d] * mask * hit_mask))
                    runs.append(run)

# plot
fig, ax = plt.subplots(1, 1)
fig.set_size_inches(10, 5)
fig.set_layout_engine('compressed')

ax.scatter(runs, photons, alpha=0.5, s=0.8)
ax.spines[['right', 'top']].set_visible(False)
ax.set_xlabel('run number')
ax.set_ylabel('photon counts')
ax.set_title(f"scatter plot of photon counts near optical axis for sizes: {1e9 * size_min} to {1e9 * size_max} nm")

plt.grid(visible=True, which='both', alpha = 0.3)
plt.savefig(f'{PREFIX}/scratch/log/peak_intensity_report.pdf')
#plt.show()

            
