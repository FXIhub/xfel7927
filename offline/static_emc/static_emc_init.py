import numpy as np
import h5py
from tqdm import tqdm
import pickle
import sys
import utils

np.random.seed(1)

class A():
    def __init__(self, C, L, D, I, mask, B, pixel_indices, file_index, frame_index, frame_shape, frame_slice, beta):
        self.betas = beta
        self.C = C
        self.L = L
        self.D = D
        self.I = I
        self.frame_shape = frame_shape
        self.frame_slice = frame_slice
        self.mask = mask
        
        self.most_likely_classes = []
        self.LL = []
        self.expectation_values = []
        self.iterations = 0
         
        self.pixel_indices = pixel_indices
        self.file_index = file_index
        self.frame_index = frame_index
        
        self.LR = np.empty((D, C), dtype = np.float32)
        self.P  = np.zeros((D, C), dtype = np.float32)
        #self.PT = np.empty((C, D), dtype = np.float32)
        
        self.w = np.ones((D,), dtype = np.float32)
        self.b = np.ones((D, L), dtype = np.float32)
        self.W = 1e-3 + np.ascontiguousarray(np.random.random((C, I)).astype(np.float32))
        self.B = np.zeros((L, I), dtype = np.float32)
        
        if type(B) is not type(None) :
            self.B[0] = B

def init(c):
    output = c.working_dir + '/recon.pickle'
    output_photons = c.working_dir + '/photons.pickle'

    if type(c.data) is str :
        fnams = [c.data]
    else :
        fnams = c.data
    
    # load mask    
    print('loading mask:', fnams[0], '/entry_1/instrument_1/detector_1/good_pixels')
    with h5py.File(fnams[0]) as f:
        mask0 = f['/entry_1/instrument_1/detector_1/good_pixels'][()]
        frame_shape = mask0.shape
        frame_size  = mask0.size
    
    # now set all pixels not selected in config to False
    mask = np.zeros_like(mask0)
    #mask.ravel()[c.pixels] = mask0.ravel()[c.pixels]
    mask[c.frame_slice] = mask0[c.frame_slice]
    print('**** number of unmasked pixels in static emc mask', mask.shape, mask.dtype, np.sum(mask))
    
    # store the un-masked pixel indices (flattened) 
    # which will be used for the reconstruction
    pixel_indices = np.arange(frame_size)[mask.ravel()]
    
    I = np.sum(mask)
    
    # load data into sparse array
    K    = []
    inds = []
    file_index = []
    frame_index = []
    inds_f = np.arange(I, dtype = np.int64)
    filter_count = 0
    low_count = 0
    
    for i, fnam in enumerate(fnams):
        if len(K) < c.max_frames :
            with h5py.File(fnam) as f:
                data = f['entry_1/data_1/data']
                # assume already applied
                # mask_per_pattern = f['/entry_1/instrument_1/detector_1/mask']
                D    = data.shape[0] 
                
                if c.filter_by is not None :
                    filter = f[c.filter_by][()]
                    filter_value = c.filter_value
                    print(f'using filter {c.filter_by} with filter value {c.filter_value}')
                    print(f'found {np.sum(filter == filter_value)} events')
                else :
                    filter = np.ones((D,), dtype=bool)
                    filter_value = True
                
                ds = np.where(filter == filter_value)[0]
                filter_count += D - len(ds)
                for d in tqdm(ds, desc = f'loading data from {fnam}'):
                    # hack for per-pattern mask
                    frame = data[d][mask].ravel() #* mask_per_pattern[d][mask].ravel()
                    m = frame > 0 
                    if (np.sum(frame[m])) > 10 :
                        K.append(frame[m].copy())
                        inds.append(inds_f[m].copy())
                        file_index.append(i)
                        frame_index.append(d)
                        
                        if len(K) == c.max_frames :
                            break
                    else :
                        low_count += 1
                        print('Warning. Frame', d, 'in dataset', fnam, 'has fewer than 10 photon counts in selected pixels')
    
    # load background 
    print('loading background:', fnams[0], '/entry_1/instrument_1/detector_1/background')
    with h5py.File(fnams[0]) as f:
        B = f['/entry_1/instrument_1/detector_1/background'][()][mask].ravel()
    
    D = len(K)
    print(f'Found {D} frames with {I} unmasked pixels')
    print(f'{low_count} frames were rejected because they had fewer than 10 photon counts in the unmasked pixels')
    print(f'{filter_count} frames were rejected due to the user defined filter dataset in the cxi files {c.filter_by}')
            
    a = A(c.classes, c.background_classes, D, I, mask, B, pixel_indices, file_index, frame_index, frame_shape, c.frame_slice, c.betas[0],)
    
    # save sparse datasets        
    print('saving reconstruction variables to:', output)
    pickle.dump(a, open(output, 'wb'))
    
    print('saving sparse photon counts to:', output_photons)
    pickle.dump([K, inds], open(output_photons, 'wb'))
    

if __name__ == '__main__' :
    config = utils.load_config(sys.argv[1])
    init(config)
