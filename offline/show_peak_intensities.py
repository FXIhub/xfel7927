import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
import glob
from tqdm import tqdm

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


runs = {}
photons = {}

for fnam in tqdm(fnams):
    run = int(fnam[-13:-9])
    with h5py.File(fnam) as f:
        if '/entry_1/sizing' in f:
            sizes  = f['/entry_1/sizing/long_axis_diameter'][()] 
            sizes += f['/entry_1/sizing/short_axis_diameter'][()] 
            sizes /= 2
               
            mask  = f['/entry_1/instrument_1/detector_1/good_pixels'][()] 
            pulse_energy = f['/entry_1/instrument_1/source_1/pulse_energy'][:]
            energy_ave = np.mean(pulse_energy)
            data = f['/entry_1/instrument_1/detector_1/data']
            
            ds = np.where( (sizes > size_min) * (sizes < size_max) )[0]
            
            if len(ds) > 0 :
                name = f['entry_1/sample_1/name'][()].decode('utf-8')
                 
                key = '/entry_1/instrument_1/detector_1/hit_score'
                if key not in f:
                    key = '/entry_1/instrument_1/detector_1/lit_pixels'
                
                lit = f[key][()][ds]

                if name not in runs :
                    runs[name] = []
                
                if name not in photons :
                    photons[name] = []
                 
                for d in ds[np.argsort(lit[::-1])][:100] :
                    photons[name].append(np.sum(data[d] * mask * hit_mask)/energy_ave)
                    runs[name].append(run)

# plot
fig, ax = plt.subplots(1, 1)
fig.set_size_inches(10, 5)
fig.set_tight_layout(True)

for name in runs.keys():
    ax.scatter(runs[name], photons[name], alpha=0.3, s=1.0, label=name)
ax.legend()
ax.set_yscale('log')
ax.spines[['right', 'top']].set_visible(False)
ax.set_xlabel('run number')
ax.set_ylabel('photon counts')
ax.set_title(f"scatter plot of photon counts near optical axis for sizes: {1e9 * size_min} to {1e9 * size_max} nm")

plt.grid(visible=True, which='both', alpha = 0.3)
plt.savefig(f'{PREFIX}/scratch/log/peak_intensity_report.pdf')
#plt.show()

            
