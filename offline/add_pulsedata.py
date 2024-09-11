# thanks to Kartik
import os
import argparse
import glob

import numpy as np
import h5py
import utils

from constants import PREFIX
from constants import XGM_DATASET, XGM_DA_NUM, WAV_DATASET, WAV_DA_NUM, TRAINID_FALLBACK

def set_values(run, da_num, dset_name, outf, out_dset_name, force=False):
    if out_dset_name in outf:
        if force:
            del outf[out_dset_name]
        else:
            print('Skipping', out_dset_name)
            return

    flist = sorted(glob.glob(PREFIX+'/raw/r%.4d/*DA%.2d*.h5'%(run, da_num)))

    with h5py.File(PREFIX + '/scratch/vds/r%.4d.cxi'%run, 'r') as f:
        vds_tid = f['entry_1/trainId'][:]
        vds_cid = f['entry_1/cellId'][:,0]
    
    out_dset = np.zeros(shape=vds_tid.shape, dtype=float)
    da_tid_dset_name = os.path.dirname(dset_name) + '/trainId'

    num_pos = 0
    for fname in flist:
        with h5py.File(fname, 'r') as f:
            if da_tid_dset_name in f :
                da_tid = f[da_tid_dset_name][:]
            else :
                da_tid = f[TRAINID_FALLBACK][:]
            # AGIPD hack for empty cell ID 0
            #vals = np.insert(f[dset_name][:], 0, 0, axis=1)

            # get unit
            prefix = f[dset_name].attrs['metricPrefixSymbol']
            
            # convert to SI
            vals = f[dset_name][:] * utils.prefix[prefix]
        
        for i, tid in enumerate(da_tid):
            pos = np.where(vds_tid == tid)[0]
            
            # if vals is pulse resolved
            if vals.ndim == 2 :
                out_dset[pos] = vals[i, vds_cid[pos]]
            
            # if vals is train resolved
            elif vals.ndim == 1 :
                out_dset[pos] = vals[i]
            
            num_pos += len(pos)
    
    utils.update_h5(outf, out_dset_name, out_dset, compression=True)

    if vds_tid.shape[0] != num_pos:
        print('WARNING: Unfilled %s values: %d vs %d' % (out_dset_name, vds_tid.shape[0], num_pos))

def main():
    parser = argparse.ArgumentParser(description='Add pulse-resolved metadata to events file')
    parser.add_argument('run', help='Run number', type=int)
    parser.add_argument('-f', '--force', help='Replace existing data if exists', action='store_true')
    args = parser.parse_args()
    
    outf = h5py.File(PREFIX + '/scratch/events/r%.4d_events.h5' % args.run, 'a')
    
    set_values(args.run, XGM_DA_NUM, XGM_DATASET, outf, 'pulse_energy', force=args.force)
    
    # probably not train resulved in reality, need to see if we can get pulse resolved wavelength
    set_values(args.run, WAV_DA_NUM, WAV_DATASET, outf, 'wavelength', force=args.force)
    
    outf.close()

if __name__ == '__main__':
    main()
