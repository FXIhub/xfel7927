import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
import glob
from tqdm import tqdm

PREFIX = os.environ["EXP_PREFIX"]

sample = 'cube'
size_min = 20e-9
size_max = 30e-9
subtract_background = False
normalise_pulse_energy = False
display_background = True


out = f'{PREFIX}/scratch/log/peak_intensity_report.pdf'
out_pickle = f'{PREFIX}/scratch/log/peak_intensity_report.pickle'

# load hit_finding mask
fnam = f'{PREFIX}/scratch/det/hit_finding_mask.h5'
with h5py.File(fnam) as f:
    hit_mask = f['entry_1/good_pixels'][()]

# loop over runs in saved hits
fnams = sorted(glob.glob(f'{PREFIX}/scratch/saved_hits/*.cxi'))


runs = {}
photons = {}
back_line = {}
average_pulse_energy = {}

back_key = '/entry_1/instrument_1/detector_1/background'
pulse_key = '/entry_1/instrument_1/source_1/pulse_energy'

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
            

            if normalise_pulse_energy and pulse_key in f:
                pulse = f[pulse_key][()]
                ds = np.where( (sizes > size_min) * (sizes < size_max) * (pulse>1e-3))[0]
                average_pulse_energy[run] = np.mean(pulse[pulse > 1e-3])
            else :
                pulse = None
                ds = np.where( (sizes > size_min) * (sizes < size_max) )[0]

            if back_key in f :
                b = f[back_key][()]
                back = np.sum(b * mask * hit_mask)
            else :
                back = None
            
            if subtract_background :
                if back_key not in f :
                    continue

            
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
                 
                back_line[run] = back
                
                for d in ds[np.argsort(lit[::-1])][:100] :
                    photons[name].append(np.sum(data[d] * mask * hit_mask)/energy_ave)
                    runs[name].append(run)

                    if subtract_background :
                        photons[name][-1] -= back
                    
                    if normalise_pulse_energy and pulse is not None and pulse[d] > 0 :
                        photons[name][-1] /= (3e3 * pulse[d])


# save
import pickle
pickle.dump((runs, photons), open(out_pickle, 'wb'))

# plot
fig, ax = plt.subplots(1, 1)
fig.set_size_inches(10, 5)
fig.set_tight_layout(True)

for name in runs.keys():
    ax.scatter(runs[name], photons[name], alpha=0.6, s=1.0, label=name)

if display_background :
    r = back_line.keys()
    if normalise_pulse_energy :
        v = [back_line[i]/(3e3 * average_pulse_energy[i]) if back_line[i] is not None else None for i in r]
    else :
        v = [back_line[i] for i in r]
    ax.scatter(r, v, alpha=0.6, c = 'k', s=1.0, label='background')

ax.legend(markerscale=5)
ax.set_yscale('log')
ax.spines[['right', 'top']].set_visible(False)
ax.set_xlabel('run number')
if subtract_background :
    ax.set_ylabel('photon counts above background')
elif normalise_pulse_energy :
    ax.set_ylabel('photon counts / 3 mJ')
elif normalise_pulse_energy and subtract_background :
    ax.set_ylabel('photon counts above background / mJ')
else :
    ax.set_ylabel('photon counts')
ax.set_title(f"scatter plot of photon counts near optical axis for sizes: {int(1e9 * size_min)} to {int(1e9 * size_max)} nm\nbackground subtraction = {subtract_background}\npulse energy normalisation = {normalise_pulse_energy}", fontsize=12)

plt.grid(visible=True, which='both', alpha = 0.3)
plt.savefig(out)
#plt.show()

            
