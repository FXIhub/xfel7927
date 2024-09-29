# fix names in cxi files due to spaces truncating name
# read run table 
# replace /entry_1/sample_1/name with name in run_table[run]['Sample']


import h5py
import glob
import os 
import json

PREFIX = os.environ['EXP_PREFIX']

# load run table
fnam = f'{PREFIX}/scratch/log/run_table.json'
with open(fnam, 'r') as f:
    run_table = json.load(f)

skey = '/entry_1/sample_1/name'
for run in run_table.keys():
    cxi_file = f'{PREFIX}/scratch/saved_hits/r{run:>04}_hits.cxi'
    print(run, cxi_file)
    if os.path.exists(cxi_file) :
        # read write file must exist
        with h5py.File(cxi_file, 'r+') as f:
            if skey in f:
                out = run_table[run]['Sample'].replace(' ', '_')
                print('replacing', f[skey][...], 'with', out, 'in file', cxi_file)
                f[skey][...] = out

    
