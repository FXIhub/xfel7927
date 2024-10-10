"""
set median as the centre of normal distribution

estimate standard deviation from left side of distribution

std = 1.4826 x (median - q25)

q25 = np.percentile(h, 25)

see: 
https://en.wikipedia.org/wiki/Interquartile_range
https://en.wikipedia.org/wiki/Median_absolute_deviation

Then apply sigma threshold 

"""
# determine threshold
import argparse

parser = argparse.ArgumentParser(description='Load hit-scores from events.h5 file and write "is_hit" and "is_miss" datasets.')
parser.add_argument('run', type=int, nargs = '+', help='Run number/s')
parser.add_argument('-t', '--hit_score_threshold_sigma',
                    help='sigma threshold for hitscore to determine is_hit',
                    type=float, default=4)
parser.add_argument('--per_train', action='store_true', help='estimate distribution per train rather than accross the entire run')
args = parser.parse_args()


import h5py
import numpy as np
import utils
from constants import PREFIX

sigma_threshold  = args.hit_score_threshold_sigma

# perhaps we set this as a percentile
# can have a low hitscore due to He 
# can have a low hitscore due to no beam 
# can have a low hitscore due to bad cell readings
# can't think of a good way to automate
minimum_hitscore = 20

for run in args.run :
    fnam = f'{PREFIX}/scratch/events/r{run:>04}_events.h5'
    #fnam = '/home/andyofmelbourne/Documents/2024/p7927/scratch/events/r0427_events.h5'
    
    if not os.path.exists(fnam) and len(args.run) > 0 :
        print('WARNING {fnam} does not exist, skipping run {run}')
        continue

    with h5py.File(fnam) as f:
        hitscore = np.rint(f['/total_intens'][()]).astype(int)
        trainIds = f['trainId'][()]
    
    unique_trainIds = np.unique(trainIds)
    
    N = hitscore.shape[0]
    
    hit_sig = np.zeros((N,), dtype = float)
    is_hit  = np.zeros((N,), dtype = bool)
    is_miss = np.zeros((N,), dtype = bool)
    
    median     = np.zeros((N,), dtype = int)
    threshold  = np.zeros((N,), dtype = float)
    threshold[:] = None
    
    # slow for now (probably doesn't matter)
    for train in unique_trainIds:
        m = np.where(train == trainIds)[0]
        
        h   = hitscore[m]
        med = np.median(h)
        std = 1.4826 * (med - np.percentile(h, 25))
        hits = h > (med + sigma_threshold * std)
        
        if med > minimum_hitscore :
            is_hit[m]    = hits
            hit_sig[m]   = (h - med)/std
            threshold[m] = (med + sigma_threshold * std)
        
        median[m]    = med
    
    hits = np.sum(is_hit)
    
    # randomly choose N indices from non-hits and 
    # trains where the threshold is valid
    miss_inds = np.random.choice(np.where(~is_hit * ~np.isnan(threshold))[0], hits, replace=False)
    is_miss[miss_inds] = True
    
    misses = np.sum(is_miss)
    
    print(f'run    : {run}')
    print(f'hits   : {hits :>10}{hits/N:>10.2%}')
    print(f'misses : {misses :>10}{misses/N:>10.2%}')
    print(f'frames : {N :>10}')
    print()
        
    out = {'is_hit': is_hit, 'is_miss': is_miss, 'hit_sigma': hit_sig}
    
    print('writing to', fnam)
    with h5py.File(fnam, 'a') as f:
        for k, v in out.items():
            utils.update_h5(f, k, v, compression=True)

"""
import pyqtgraph as pg
# plot result for one train
hist = np.bincount(h)
bins = np.arange(hist.shape[0])
plot = pg.plot(bins, hist)
plot.plot(bins, hist.max()/2 * np.exp( - (med - bins)**2 / (2 * std**2)), pen=pg.mkPen('r'))

plot = pg.plot(hitscore)
plot.plot(median, pen=pg.mkPen('y'))
plot.plot(threshold, pen=pg.mkPen('g'))
    
hist = np.bincount(hitscore)
plot = pg.plot(hist)
"""
