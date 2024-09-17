# determine threshold
import argparse

parser = argparse.ArgumentParser(description='Load hit-scores from events.h5 file and write "is_hit" and "is_miss" datasets.')
parser.add_argument('run', type=int, help='Run number')
parser.add_argument('-t', '--hit_score_threshold_sigma',
                    help='sigma threshold for hitscore to determine is_hit',
                    type=float, default=3)
#parser.add_argument('-m', '--modules',
#                    nargs='+',
#                    help='modules to include in hitscore calculation',
#                    type=int, default=[4,])
args = parser.parse_args()

import utils
from constants import PREFIX
from scipy import optimize
import numpy as np
import h5py

events_fnam = f'{PREFIX}/scratch/events/r{args.run:>04}_events.h5'

# get hitscores
# use integrated photon counts
with h5py.File(events_fnam) as f:
    hitscore = np.sum(f['total_intens'][:], axis=1)

# Get hit indices
def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x-x0)**2 / 2 / sigma**2)

hit_score_threshold_sigma = 3

offset = 100
hy, hx = np.histogram(hitscore, np.arange(0, hitscore.max()+1, 1))
xmax   = hy[offset:].argmax() + offset # Ignoring first few bins
hymax  = hy[xmax]

m = (hitscore > offset) * (hitscore < (2*xmax - offset))
sig = np.std(hitscore[m])

popt, pcov = optimize.curve_fit(gaussian, hx[offset:xmax], hy[offset:xmax], p0=(hy.max(), xmax, sig))

# thresholds for miss
miss_thresh_max = popt[1] - np.abs(args.hit_score_threshold_sigma   * popt[2])
miss_thresh_min = popt[1] - 3*np.abs(args.hit_score_threshold_sigma * popt[2])

# thresholds for hit
hit_thresh_min = popt[1] + np.abs(hit_score_threshold_sigma * popt[2])
hit_thresh_max = np.inf

print('Fitted background Gaussian to hit score: %.3f +- %.3f' % (popt[1], popt[2]))
print('Applying a hitscore threshold of       : %.3f' % (hit_thresh_min))

# write to events file
is_hit  = (hitscore > hit_thresh_min)  * (hitscore < hit_thresh_max)
is_miss = (hitscore > miss_thresh_min) * (hitscore < miss_thresh_max)

print('\nFound a total of:')
N = len(hitscore)
hits = np.sum(is_hit)
misses = np.sum(is_miss)
print(f'hits   : {hits :>10}{hits/N:>10.2%}')
print(f'misses : {misses :>10}{misses/N:>10.2%}')
print(f'frames : {N :>10}')
print()

out = {'is_hit': is_hit, 'is_miss': is_miss}

print('writing to', events_fnam)
with h5py.File(events_fnam, 'a') as f:
    utils.update_h5(f, 'is_hit', is_hit, compression=True)
    utils.update_h5(f, 'is_miss', is_miss, compression=True)

# plot 
#import pyqtgraph as pg
#plot = pg.plot(hx[1:], hy)
#plot.plot(hx, gaussian(hx, popt[0], popt[1], popt[2]))

