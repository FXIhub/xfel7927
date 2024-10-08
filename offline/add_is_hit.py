# determine threshold
import argparse

parser = argparse.ArgumentParser(description='Load hit-scores from events.h5 file and write "is_hit" and "is_miss" datasets.')
parser.add_argument('run', type=int, nargs = '+', help='Run number/s')
parser.add_argument('-t', '--hit_score_threshold_sigma',
                    help='sigma threshold for hitscore to determine is_hit',
                    type=float, default=3)
parser.add_argument('--per_train', action='store_true', help='estimate distribution per train rather than accross the entire run')
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
from tqdm import tqdm

for run in args.run :
    events_fnam = f'{PREFIX}/scratch/events/r{run:>04}_events.h5'
    
    # Get hit indices
    def gaussian(x, a, x0, sigma):
        return a * np.exp(-(x-x0)**2 / 2 / sigma**2)
    
    def get_hits_misses(hitscore, offset = 10):
        hy, hx = np.histogram(hitscore, np.arange(0, hitscore.max()+1, 1))
        xmax   = hy[offset:].argmax() + offset # Ignoring first few bins
        hymax  = hy[xmax]
        
        m = (hitscore > offset) * (hitscore < (2*xmax - offset))
        sig = np.std(hitscore[m])
        
        popt, pcov = optimize.curve_fit(gaussian, hx[offset:xmax], hy[offset:xmax], p0=(hy[offset:].max(), xmax, sig))

        mean, sig = popt[1], popt[2]

        # thresholds for miss -3 to +3 sigma from mean
        t = args.hit_score_threshold_sigma - 1
        miss_thresh_max = mean + t*np.abs(sig)
        miss_thresh_min = mean - t*np.abs(sig)
        
        # thresholds for hit
        hit_thresh_min = mean + np.abs(args.hit_score_threshold_sigma * sig)
        hit_thresh_max = np.inf
        
        #print('Fitted background Gaussian to hit score: %.3f +- %.3f' % (popt[1], popt[2]))
        #print('Applying a hitscore threshold of       : %.3f' % (hit_thresh_min))
        
        # write to events file
        is_hit  = (hitscore > hit_thresh_min)  * (hitscore < hit_thresh_max)
        is_miss = (hitscore > miss_thresh_min) * (hitscore < miss_thresh_max)
        
        return is_hit, is_miss
    
    def get_hits_misses_simple(hitscore, offset=10):
        m = (hitscore > offset) * (hitscore < hitscore.max())
           
        mean, sig = np.mean(hitscore[m]), np.std(hitscore[m])
        
        # thresholds for miss -3 to 0 sigma from mean
        miss_thresh_max = mean 
        miss_thresh_min = mean - 3*np.abs(sig)
        
        # thresholds for hit
        hit_thresh_min = mean + np.abs(args.hit_score_threshold_sigma * sig)
        hit_thresh_max = np.inf
        
        #print('Fitted background Gaussian to hit score: %.3f +- %.3f' % (popt[1], popt[2]))
        #print('Applying a hitscore threshold of       : %.3f' % (hit_thresh_min))
        
        # write to events file
        is_hit  = (hitscore > hit_thresh_min)  * (hitscore < hit_thresh_max)
        is_miss = (hitscore > miss_thresh_min) * (hitscore < miss_thresh_max)
        
        return is_hit, is_miss

    # get hitscores
    # use integrated photon counts
    with h5py.File(events_fnam) as f:
        hitscore = f['litpixels'][:]

    if not args.per_train :
        is_hit, is_miss = get_hits_misses(hitscore)
    else :
        # determine per N-events (1s)
        N = 10 * 352 
        
        is_hit  = np.zeros(hitscore.shape[0], dtype = bool)
        is_miss = np.zeros(hitscore.shape[0], dtype = bool)
        index = 0
        for i in tqdm(range(0, hitscore.shape[0], N), disable = False):
            start = i 
            stop  = i+N 
            
            h, m = get_hits_misses_simple(hitscore[start: stop])
            is_hit[start: stop]  = h
            is_miss[start: stop] = m
            
        

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

