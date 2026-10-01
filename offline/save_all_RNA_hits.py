
import numpy as np
import os 

PREFIX = os.environ['EXP_PREFIX']

cxi_out = f'{PREFIX}/scratch/saved_hits/RNA_all_hits.cxi'

runs = [619,620,621,622,623,624,625,626,627,628,629]

cxi_in  = []
for run in runs :
    cxi_in.append(f'{PREFIX}/scratch/saved_hits/r{run:>04}_hits.cxi')

def cxi_filter(f):
    out = np.ones((f['/entry_1/data_1/data'].shape[0],), dtype = bool)
    return out


