from pathlib import Path
import numpy as np
import glob
import argparse
import h5py
import utils
import os
import json

PREFIX = os.environ["EXP_PREFIX"]

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter, 
    description="""
Load the run table and return a list of runs with a given attribute/s. 

For example to get a list of all gold runs:
    python get_runs.py -k Sample -v AuNP
    >>> 1,2,3

For mulitple conditions append keys (-k) and values (-k) in order.

So to get all Ery runs where the cxi files have been processed:
    python get_runs.py -k Sample CXI -v Erythrocruorin ready
    >>> 51,52,55

Use the spreadsheet https://docs.google.com/spreadsheets/d/15nfDVdlXoTenMnGrbUlyC005MD8xmpgZX0R5y1fnQeI/edit?gid=85516743#gid=85516743
to look for key and Value pairs

This is useful for submiting slurm jobs:
    submit_events.sh $(python get_runs.py -k Sample -v Erythrocruorin)
""")

default_run_table = f'{PREFIX}scratch/log/run_table.json'

parser.add_argument('-k', '--keys',   type=str, nargs='+', required = True, help='key')
parser.add_argument('-v', '--values', type=str, nargs='+', required = True, help='value')
parser.add_argument('--no_dashes', action='store_true', help='do not use dashes to join consecutive run numbers together e.g. 1,3,6,7,8,9,15 -> 1,2,6-9,15')
parser.add_argument('--join', type=str, default = ',', help='str used to join numbers together')
parser.add_argument('--format', type=str, help="replace run numbers with another string e.g.: --format 'r{run:>04}_powder.h5' produces r0001_powder.h5,r0004_powder.h5. This disables dashes.")
parser.add_argument('-r', '--run_table', 
                    help='Location of the run_table json file',
                    default=default_run_table)
args = parser.parse_args()

# newlines are escaped by argparse
# but we want to allow them
args.join = args.join.replace('\\n', '\n')

assert(len(args.keys) == len(args.values))

# load run table
run_table = json.load(open(args.run_table, 'r'))

runs = []
for run, run_dict in run_table.items():
    if not isinstance(run_dict, dict):
        continue
    
    add = True
    for k, v in zip(args.keys, args.values) :
        if k not in run_dict:
            add = False
            break

        if v not in run_dict[k]:
            add = False
            break

    if add :
        runs.append(run_dict['Run number'])

runs = sorted(runs)

# use - for consecutive number ranges
out = []
if not args.no_dashes and not args.format:
    # if this is the start of sequence > 2 of consecutive numbers
    i = 0
    while i < len(runs):
        # test for begining of sequence
        j = i + 1 
        while j < len(runs) :
            if runs[j] == (runs[j-1] + 1) :
                j += 1
            else :
                break
        if (j-i) <= 2 :
            out.append(f'{runs[i]}')
            i += 1
        else :
            out.append(f'{runs[i]}-{runs[j-1]}')
            i = j
elif args.format :
    out = [args.format.format(run=run) for run in runs]
else :
    out = [str(run) for run in runs]


a = args.join.join(out)
print(a)
