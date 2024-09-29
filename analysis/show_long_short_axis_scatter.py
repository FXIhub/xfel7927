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

#out = f'{PREFIX}/scratch/log/peak_intensity_report.pdf'
#out_pickle = f'{PREFIX}/scratch/log/peak_intensity_report.pickle'

# load hit_finding mask
fnam = f'{PREFIX}/scratch/det/hit_finding_mask.h5'
with h5py.File(fnam) as f:
    hit_mask = f['entry_1/good_pixels'][()]

# loop over runs in saved hits
fnams = sorted(glob.glob(f'{PREFIX}/scratch/saved_hits/*.cxi'))

long  = {}
short = {}

for fnam in tqdm(fnams):
    run = int(fnam[-13:-9])
    with h5py.File(fnam) as f:
        if '/entry_1/sizing' in f:
            name       = f['entry_1/sample_1/name'][()].decode('utf-8')
            
            if name not in long:
                long[name] = []
            if name not in short:
                short[name] = []
            
            l = f['/entry_1/sizing/long_axis_diameter'][()] 
            s = f['/entry_1/sizing/short_axis_diameter'][()] 
            long[name] += [i for i in l]
            short[name] += [i for i in s]
                 
# save
#import pickle
#pickle.dump((runs, photons), open(out_pickle, 'wb'))

# plot
fig, ax = plt.subplots(1, 1)
fig.set_size_inches(10, 5)
fig.set_tight_layout(True)

for name in long.keys():
    if 'Au' not in name:
        ax.scatter(1e9 * np.array(long[name]), 1e9 * np.array(short[name]), alpha=0.7, s=4.0, label=name)

ax.legend(markerscale=5, fontsize=20)
#ax.set_yscale('log')
ax.spines[['right', 'top']].set_visible(False)
ax.set_xlabel('long axis diameter (nm)')
ax.set_ylabel('short axis diameter (nm)')
ax.set_title(f"scatter plot of hit sizes")

#plt.grid(visible=True, which='both', alpha = 0.3)
#plt.savefig(out)
plt.show()

            
