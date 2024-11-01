from pathlib import Path
import numpy as np
import glob
import argparse
import h5py
import utils
import os

PREFIX = os.environ["EXP_PREFIX"]

def main():
    parser = argparse.ArgumentParser(description='calculate powder per pixel for all per-cell powder patterns matching r<run>_powder_*.h5')
    parser.add_argument('run', type=int, help='Run number')
    parser.add_argument('-o', '--powder_folder', 
                        help='Path of powder folder (default=%s/scratch/powder/)'%PREFIX,
                        default=PREFIX+'/scratch/powder/')
    args = parser.parse_args()
    
    mask_fnam = f'{PREFIX}/scratch/det/r{args.run:>04}_mask.h5'
    
    # get matching per cell powder patterns
    fnams = glob.glob(f'{args.powder_folder}r{args.run:>04}_powder_*.h5')
    fnams = [fnam for fnam in fnams if 'pixel' not in fnam]
    
    # out is fnam + 'per_pixel'
    out_fnams = []
    for fnam in fnams :
        p = Path(fnam)
        out = p.with_stem(p.stem + '_per_pixel')
        out_fnams.append(str(out))

    # load mask
    print(f'loading per-cell mask {mask_fnam}')
    with h5py.File(mask_fnam) as f:
        cellIds = f['entry_1/cellIds'][()]
        mask = {}
        for i, c in enumerate(cellIds):
            mask[c] = f['entry_1/good_pixels'][i]
    
    for out, fnam in zip(out_fnams, fnams) :
        print(f'processing powder {fnam}')
        with h5py.File(fnam) as f:
            modules      = f['modules'].shape[0]
            module_shape = f['data'].shape[-2:]
            shape = (modules,) + f['data'].shape[-2:]
        
        powder      = np.zeros(shape, dtype = float)
        powder_sum  = np.zeros(module_shape, dtype = np.uint64)
        overlap     = np.zeros(module_shape, dtype = np.uint64)
        
        for module in range(modules):
            # load powder for module for all cells
            with h5py.File(fnam) as f:
                powder_cell = f['data'][module]
                events_cell = f['events'][module]
                cellIds     = f['cellIds'][()]
            
            powder_sum.fill(0)
            overlap.fill(0)
            
            for i, c in enumerate(cellIds):
                powder_sum += (mask[c][module] * powder_cell[i]).astype(powder_sum.dtype)
                overlap    += (mask[c][module] * events_cell[i]).astype(overlap.dtype)
            
            powder[module] = powder_sum / np.clip(overlap, 1, None)
            
        print(f'writing per pixel powder to: {out}')
        with h5py.File(out, 'w') as f:
            utils.update_h5(f, 'data', powder.astype(np.float32), compression = True, chunks = powder.shape)


if __name__ == '__main__':
    main()
