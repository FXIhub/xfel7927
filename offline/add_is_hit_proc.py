# determine threshold
import argparse

parser = argparse.ArgumentParser(description='Use extra_redu to perform hitfinding and write "is_hit" and "is_miss" datasets to events file. need: pip install extra-redu')
parser.add_argument('run', type=int, nargs = '+', help='Run number/s')
parser.add_argument('-t', '--hit_score_threshold_sigma',
                    help='sigma threshold for hitscore to determine is_hit',
                    type=float, default=4)
parser.add_argument('-m', '--modules',
                    nargs='+',
                    help='modules to include in hitscore calculation',
                    type=int, default=[15])
args = parser.parse_args()



# read litpixels object from extra data 
# change hitfinding parameters
# hit find 
# then write to events file
import numpy as np
import h5py

from extra_redu.fileutils import StackedPulseSource
from extra_data import open_run
from extra_redu.spi import SpiHitfinder

import utils
from constants import PREFIX, EXP_ID

for run in args.run :
    events_fnam = f'{PREFIX}/scratch/events/r{run:>04}_events.h5'
    run = open_run(EXP_ID, run, data='all')
    
    # hitfindig results
    run['SPB_DET_AGIPD1M-1/REDU/SPI_HITFINDER:output'].keys()
    # data.hitFlag - boolean indicates hits
    # data.missFlag - boolean indicates selected misses for background characterization
    # data.dataFlag = data.hitFlag | data.missFlag
    # data.hitscore - hitscore
    # data.pulseId - pulse Id
    # data.trainId - train Id

    # hitfinder configuration
    #a = run['SPB_DET_AGIPD1M-1/REDU/SPI_HITFINDER']
    #for k in a.keys():
    #    print(k, a[k].ndarray()[0])

    # litpixels
    #a = run['SPB_DET_AGIPD1M-1/CORR/0CH0:output'].keys()
    #print(a)
    #for k in a.keys():
    #    print(k, a[k].ndarray()[0])

    #src = StackedPulseSource.from_datacollection(
    #    run, r"SPB_DET_AGIPD1M-1/DET/(?P<key>\d+)CH0:xtdf", "image")

    litpx_src = StackedPulseSource.from_datacollection(
            run, f"SPB_DET_AGIPD1M-1/CORR/(?P<key>\d+)CH0:output", "litpx")

    # module 3 has a bad panel
    hitfinder = SpiHitfinder(
        modules=args.modules,
        mode="adaptive",
        snr=args.hit_score_threshold_sigma,
        min_scores=100,
        fixed_threshold=0,
        hitrate_window_size=200,
        miss_fraction=1,
        miss_fraction_base="hit",
        xgm_norm = False,
    )
    hitfinder.find_hits(litpx_src)

    hit_score = hitfinder.hitscore
    is_hit  = hitfinder.hit_mask
    is_miss = hitfinder.miss_mask
    
    print('\nFound a total of:')
    N = len(is_hit)
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
        utils.update_h5(f, 'hit_score', hit_score, compression=True)



