import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import os
import glob
from tqdm import tqdm


import sys
run_non_interactive = 'background' in sys.argv[1]


PREFIX = os.environ["EXP_PREFIX"]

sample = 'cube'
size_min = 20e-9
size_max = 30e-9
subtract_background = False
normalise_pulse_energy = False
display_background = True

# just before refinement
# + 0.5 so that we do not overlap with scatter
focus_refined = np.array([111, 197, 285, 333, 371, 404, 413, 436, 456, 468, 494, 500, 536]) + 0.5

# where He was not used inclusive
no_he = [[1, 72], [111, 124], [221, 234], [336, 348], [484, 495]]


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
files = {}
indexes = {}

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
                print(fnam)
                name = f['entry_1/sample_1/name'][()].decode('utf-8')
                 
                key = '/entry_1/instrument_1/detector_1/hit_score'
                if key not in f:
                    key = '/entry_1/instrument_1/detector_1/lit_pixels'
                
                lit = f[key][()][ds]
                
                if name not in runs :
                    runs[name] = []
                
                if name not in photons :
                    photons[name] = []
                
                if name not in files :
                    files[name] = []
                    indexes[name] = []
                 
                back_line[run] = back
                
                for d in ds[np.argsort(lit[::-1])][:100] :
                    photons[name].append(np.sum(data[d] * mask * hit_mask))
                    runs[name].append(run)
                    files[name].append(fnam)
                    indexes[name].append(d)
                    
                    if subtract_background :
                        photons[name][-1] -= back
                    
                    if normalise_pulse_energy and pulse is not None and pulse[d] > 0 :
                        photons[name][-1] /= (3e3 * pulse[d])


def load_frame(fnam, index):
    with h5py.File(fnam) as f:
        frame = f['entry_1/data_1/data'][index]
        #image = geom.position_modules(frame)[0]
    return frame


if run_non_interactive == False :
    #import pyqtgraph as pg
    fig_im, ax_im = plt.subplots()
    fig_im.set_size_inches(10, 10)
    fig_im.set_tight_layout(True)

    import extra_geom
    geom_fnam = '../geom/r0300.geom'
    geom = extra_geom.AGIPD_1MGeometry.from_crystfel_geom(geom_fnam)

    name = list(files.keys())[0]
    fnam = files[name][0]
    i    = indexes[name][0]

    #geom.plot_data(load_frame(fnam, i), ax = ax_im, colorbar=False)
    frame_plot = ax_im.imshow(geom.position_modules(load_frame(fnam, i))[0]**0.2, vmin =0)
    fig_im.show()

    def on_pick(event):
        i = event.ind[0]
        run   = runs[event.artist.name][i]
        fnam  = files[event.artist.name][i]
        index = indexes[event.artist.name][i]
        print(f'clicked on run {run}, fnam {fnam}, index {index}')
        
        frame = load_frame(fnam, index)
        
        #geom.plot_data(frame, ax = ax_im, vmax = 3)
        frame_plot.set_data(geom.position_modules(frame)[0]**0.2)
        fig_im.canvas.draw()
        fig_im.canvas.flush_events()
        #ax_im.draw()
    


# save
import pickle
pickle.dump((runs, photons), open(out_pickle, 'wb'))

# plot
fig, ax = plt.subplots(1, 1)
fig.set_size_inches(30, 8)
fig.set_tight_layout(True)

artists = []
for name in runs.keys():
    artists.append(ax.scatter(runs[name], photons[name], alpha=0.6, s=3.0, picker=5, label=name))
    artists[-1].name = name

if display_background :
    r = back_line.keys()
    if normalise_pulse_energy :
        v = [back_line[i]/(3e3 * average_pulse_energy[i]) if back_line[i] is not None else None for i in r]
    else :
        v = [back_line[i] for i in r]
    ax.scatter(r, v, alpha=0.6, c = 'k', s=3.0, label='background')

# discard low signal
ylim = ax.get_ylim()
ax.set_ylim([50, ylim[1]])

# add vlines when beam was refined
ylim = ax.get_ylim()
ax.vlines(focus_refined, ylim[0], ylim[1], label='focus refinement', color='k', linestyle='--')


# colour areas where no He was used
label = 'no He'
for i, j in no_he :
    ax.fill_between([i, j], 0, 1, 
                    color='grey', alpha=0.2, transform=ax.get_xaxis_transform(), label=label)
    label = None


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

# For the minor ticks, use no labels; default NullFormatter.
ax.xaxis.set_minor_locator(MultipleLocator(10))

if run_non_interactive == False :
    fig.canvas.callbacks.connect('pick_event', on_pick)

ax.grid(visible=True, which='both', alpha = 0.3)
fig.show()

plt.savefig(out)

if run_non_interactive == False :
    plt.show()
            

