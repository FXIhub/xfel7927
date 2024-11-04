
import argparse
import numpy as np
import h5py
from tqdm import tqdm
import os
import utils
import common
import time

from constants import PREFIX, EXP_ID, NMODULES, MODULE_SHAPE

def get_noise_cal(propno, runno, cellIds, module):
    from extra_data import open_run
    from extra_geom import AGIPD_1MGeometry
    from datetime import datetime
    from extra.calibration import CalibrationData, AGIPDConditions
    
    run = open_run(propno, runno, data="raw")

    run_date = datetime.strptime(run.run_metadata()["creationDate"], "%Y%m%dT%H%M%SZ")
    ncell = int(run["SPB_IRU_AGIPD1M1/MDL/FPGA_COMP"].run_value("bunchStructure.nPulses"))
    acquisition_rate = float(run["SPB_IRU_AGIPD1M1/MDL/FPGA_COMP"].run_value("bunchStructure.repetitionRate"))
    source_energy = 9.2
    bias_voltage = 300.0
    integration_time = float(run["SPB_IRU_AGIPD1M1/MDL/FPGA_COMP"].run_value("integrationTime"))
    gain_mode = run["SPB_IRU_AGIPD1M1/MDL/FPGA_COMP"].run_value("gainModeIndex")
    gain_setting = int(run["SPB_IRU_AGIPD1M1/MDL/FPGA_COMP"].run_value("gain"))
    
    acond = AGIPDConditions(bias_voltage, ncell, acquisition_rate, gain_setting, gain_mode, source_energy, integration_time)

    # this seems to randomly fail
    for i in range(100):
        try :
            agipd_cd = CalibrationData.from_condition(acond, 'SPB_DET_AGIPD1M-1', event_at = run_date)
            break
        except :
            time.sleep(1)
    
    # don't know if cellId's start at 0 or 1, I think 0
    # AGIPD03
    s = agipd_cd['Noise'][f'AGIPD{module:>02}'].ndarray()
    # (128, 512, 352, 3)
    gain = 0
    s = s[..., gain].transpose((2, 1, 0))
    return s[cellIds]
    
def get_noise_cal_testing(propno, runno, cellIds, module):
    with h5py.File('r0600_noise.h5', 'r') as f:
        sig = f['Noise'][cellIds, module]
    return sig

def discard_outliers(ar, sig_min, sig_max):
    """
    2 pass MAD filter
    """
    #med = np.median(ar)
    #dev = ar - med
    #MAD = np.median(np.abs(dev))
    q1 = np.percentile(ar, 25)
    q3 = np.percentile(ar, 75)
    med = (q1 + q3) / 2
    dev = ar - med
    MAD = (q3 - q1) / 2
    std = 1.4826 * MAD
    
    mask = (dev > sig_min*std) * (dev < sig_max*std)
    if np.any(mask) :
        q1 = np.percentile(ar, 25)
        q3 = np.percentile(ar, 75)
        med = (q1 + q3) / 2
        MAD = (q3 - q1) / 2
        std = 1.4826 * MAD
    
    minval = med + sig_min * std
    maxval = med + sig_max * std
    mask = (ar > minval) * (ar < maxval)
    return mask, minval, maxval, {'med': med, 'std': std}

def make_sigma_mask(sig, std_min = -3, std_max = 1.5, Nmax = 20):
    mask, minval, maxval, info = discard_outliers(sig, std_min, std_max)
    
    # if more than 20 cells have bad sigma then discard pixel
    pixel_mask = np.sum(~mask, axis=0) <= Nmax 
    mask *= pixel_mask[None, ...]
    
    print('sigma mask min, max, % masked:', minval, maxval, round(100*np.sum(~mask) / mask.size, 2), '%')
    return mask

def make_powder_mask(powder_cell, sig_max = 3):
    ar = np.ascontiguousarray(np.transpose(powder_cell, (1,2,0)))
    # discard cells with outlier counts
    q1 = np.percentile(ar, 25, axis=-1)
    q3 = np.percentile(ar, 75, axis=-1)
    med = (q1 + q3) / 2
    MAD = (q3 - q1) / 2
    std = 1.4826 * MAD
    
    maxval = med + sig_max * std
    mask   = (ar < maxval[..., None])
    
    mask = np.ascontiguousarray(np.transpose(mask, (2,0,1)))
    return mask

def make_edge_mask(shape = (512, 128)):
    edges_mask = np.ones(shape, dtype = bool)
    edges_mask[::64] = False
    edges_mask[63::64] = False
    edges_mask[:, 0] = False
    edges_mask[:, -1] = False
    return edges_mask

def masked_median_filter_2D(ar, mask, size = 8):
    s2 = size // 2
    med_fil = np.zeros_like(ar)
    for i in tqdm(range(ar.shape[0])):
        for j in range(ar.shape[1]):
            i0 = max(i-s2, 0)
            i1 = min(i+s2, ar.shape[0])
            j0 = max(j-s2, 0)
            j1 = min(j+s2, ar.shape[1])
            m = mask[i0:i1, j0:j1]
            if np.any(m) :
                med = np.median(ar[i0:i1, j0:j1][m])
            else :
                med = 0
            med_fil[i, j] = med
    return med_fil

def make_med_filter_mask(ar, mask, size = 8, ratio = 1.2):
    out = np.ones(ar.shape, dtype = bool)
    ar_fil = masked_median_filter_2D(ar, mask, size)
    m = ar_fil == 0
    ar_fil[m] = 1
    out = (ar / ar_fil) < ratio
    return out

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    Generate per-cell pixel mask based on calibrated data and per-cell powder patterns.

    First discard outlyer standard deviation values per cell in calibration 'Noise' data:
    std_sig_min, std_sig_max

    If more than 'N' cells in a pixels have outlyer std values then mask entire pixel:
    max_bad_cells

    If the counts in the per-cell powder pattern are too big compared to other cells in
    the same pixel, then mask that cell. 
    std_powder_cell

    If the per-pixel powder pattern has pixel counts much greater than it's neighbours
    then mask all cells under that pixel.
    median_filter_size
    median_filter_ratio
    """)
    parser.add_argument('run', type=int, help='Run number/s')
    parser.add_argument('--std_sig_min', type=float, default=-3)
    parser.add_argument('--std_sig_max', type=float, default=1.5)
    parser.add_argument('--max_bad_cells', type=int, default=20)
    parser.add_argument('--std_powder_cell', type=float, default=3)
    parser.add_argument('--median_filter_size', type=int, default=8)
    parser.add_argument('--median_filter_ratio', type=float, default=1.2)
    parser.add_argument('--add_global_mask', type=str, default='stray_light_mask.h5')
    args = parser.parse_args()

    args.powder = f'{PREFIX}/scratch/powder/r{args.run:>04}_powder.h5'
    args.output = f'{PREFIX}/scratch/det/r{args.run:>04}_mask.h5'
    
    if args.add_global_mask :
        args.add_global_mask = f'{PREFIX}/scratch/det/{args.add_global_mask}'
        print(f'loading user selected global mask {args.add_global_mask}')
        with h5py.File(args.add_global_mask) as f:
            global_mask = f['entry_1/good_pixels'][()]
    else :
        global_mask = np.ones((NMODULES,) + MODULE_SHAPE, dtype = bool)
    
    # get cellIds from powder
    with h5py.File(args.powder) as f:
        cellIds     = np.unique(f['cellIds'][()])
        powder_cell_dtype = f['data'].dtype
    
    # this shape helps reading per-cell 
    mask = np.ones((len(cellIds), NMODULES) + MODULE_SHAPE, dtype = bool)
    module_mask = np.ones((len(cellIds),) + MODULE_SHAPE, dtype = bool)
    
    # powder for display
    p    = np.zeros((NMODULES,) + MODULE_SHAPE, dtype = float)
    
    for module in range(NMODULES):
        print(f'processing module {module}')
        # get sigma
        #sig = get_noise_cal_testing(EXP_ID, args.run, cellIds, module)
        sig = get_noise_cal(EXP_ID, args.run, cellIds, module)
        
        # mask sigma
        module_mask = make_sigma_mask(sig, args.std_sig_min, args.std_sig_max, args.max_bad_cells)
        print('sigma mask % masked:', round(100*np.sum(~module_mask) / module_mask.size, 2), '%')
        
        # mask edges
        module_mask *= make_edge_mask(shape = MODULE_SHAPE)
        print('edge mask % masked:', round(100*np.sum(~module_mask) / module_mask.size, 2), '%')
        
        # mask powder_cells
        # data {16, 351, 512, 128}
        # cellIds {16, 351}
        print(f'reading powder from: {args.powder}')
        powder_cell = np.zeros((len(cellIds),) + MODULE_SHAPE, dtype = powder_cell_dtype)
        with h5py.File(args.powder) as f:
            cs     = f['cellIds'][module]
            for cellId in cellIds:
                c = np.where(cellId == cs)[0]
                if len(c) == 1 :
                    powder_cell[c] = f['data'][module, c]
        
        module_mask *= make_powder_mask(powder_cell, args.std_powder_cell)
        print('powder mask % masked:', round(100*np.sum(~module_mask) / module_mask.size, 2), '%')
        
        m = np.sum(module_mask, axis=0) 
        powder = np.sum(powder_cell * module_mask, axis = 0) / np.clip(m, 1, None)
        
        #module_mask *= make_med_filter_mask(powder, m>0, args.median_filter_size, args.median_filter_ratio)
        module_mask *= make_med_filter_mask(powder, m>0, args.median_filter_size, args.median_filter_ratio)
        print('median filter mask % masked:', round(100*np.sum(~module_mask) / module_mask.size, 2), '%')
        
        mask[:, module] = module_mask
         
        # testing
        #m = np.sum(module_mask, axis=0) 
        #powder = np.sum(powder_cell * module_mask, axis = 0) / np.clip(m, 1, None)
        #p[module] = powder
        #print(p.dtype)
    
    # output
    print(f'writing mask and cellids to {args.output}')
    with h5py.File(args.output, 'w') as f:
        utils.update_h5(f, 'entry_1/good_pixels', global_mask * mask, compression=True, chunks = (1,) + mask.shape[1:])
        utils.update_h5(f, 'entry_1/cellIds', cellIds, compression=True, chunks = cellIds.shape)
    
    # show powder
    #geom_fnam = common.get_geom(args.run)
    #import extra_geom
    #geom = extra_geom.AGIPD_1MGeometry.from_crystfel_geom(geom_fnam)
    #import pyqtgraph as pg
    #pg.show(geom.position_modules(p)[0])
    # allow Control-C
    #import signal
    #signal.signal(signal.SIGINT, signal.SIG_DFL) 
    #pg.exec()
