# copy a subset of patterns to data.cxi
import h5py
import numpy as np
from tqdm import tqdm
import os
import sys
import pathlib
import runpy

def load_config(path):
    p = pathlib.Path(path)
    
    # returns a dict
    config = runpy.run_path(str(p.absolute()))
    
    return config

from constants import PREFIX

#cxi_out = f'{PREFIX}/scratch/saved_hits/Ery_strong.cxi'
max_frames  = 1e10

config     = load_config(sys.argv[1])
cxi_in     = config['cxi_in']
cxi_out    = config['cxi_out']
cxi_filter = config['cxi_filter']

#cxi_in = []
# runs = range(642)
#for run in runs :
#    fnam = PREFIX+'/scratch/saved_hits/r%.4d_hits.cxi'%run
#    if os.path.exists(fnam) :
#        with h5py.File(fnam) as f:
#            name = f['/entry_1/sample_1/name'][()].decode('utf-8')
#            if 'Ery' in name :
#                cxi_in.append(fnam)


#def cxi_filter(f):
#    h = f['/entry_1/instrument_1/detector_1/hit_sigma'][()]
#    l0 = f['/entry_1/sizing/long_axis_diameter'][()]
#    s0 = f['/entry_1/sizing/short_axis_diameter'][()]
#    l  = np.where(s0 < l0, l0, s0)
#    s  = np.where(s0 < l0, s0, l0)
#
#    # long axis has to be 32nm +- error
#    # short axis (in projection) can be anywhere between 15-32nm +- error
#    out = (l >= 20e-9) * (l <= 35e-9) * (s >= 10e-9) * (s <= 35e-9) * (h > 20.)
#    #out = (h > 200.)
#    return out
        


# copy non event related data from first file
with h5py.File(cxi_in[0]) as f:
    Nevents = f['/entry_1/data_1/data'].shape[0]
    frame_shape = f['/entry_1/data_1/data'].shape[1:]
    
    with h5py.File(cxi_out, 'w') as g:
        def write_data(name, obj):
            if type(obj) == h5py.Dataset :
                print('writing', name)
                if not ((len(obj.shape) > 0) and (obj.shape[0] == Nevents)) :
                    g[name] = obj[()]
            else :
                print('skipping:', name)
        
        f.visititems(write_data)

# get total number of frames
inds   = {}
frames = 0
for fnam in cxi_in :
    with h5py.File(fnam) as f:
        inds[fnam] = np.where(cxi_filter(f))[0]
        
        frames += len(inds[fnam])

# merge good_pixels
print('merging good_pixels')
key = 'entry_1/instrument_1/detector_1/good_pixels'
mask = np.zeros(frame_shape, dtype = int)
for fnam in cxi_in :
    with h5py.File(fnam) as f:
        mask += f[key][()]
        print(f'updating mask with {round(100 * np.sum(mask==0) / mask.size, 2)}% of pixels masked', fnam)

mask = mask > 0 
with h5py.File(cxi_out, 'r+') as f:
    f[key][:] = mask
        

# stack background images 
key = '/entry_1/instrument_1/detector_1/background'
backgrounds = []
background_inds = np.empty((frames,), dtype = np.uint32)
index = 0
for i, fnam in enumerate(cxi_in) :
    with h5py.File(fnam) as f:
        Nevents = f['/entry_1/data_1/data'].shape[0]
        
        print(f'loading background from {fnam}')
        backgrounds.append(f[key][()])
        background_inds[index: index + len(inds[fnam])] = i
    
    index += len(inds[fnam])

backgrounds = np.array(backgrounds)

# write to cxi file
with h5py.File(cxi_out, 'r+') as g:
    if key in g :
        del g[key]
    g.create_dataset(key, data = backgrounds, compression = 'gzip', compression_opts = 1, chunks = (1,) + backgrounds.shape[1:])
    g.create_dataset('entry_1/background_index', data = background_inds, compression = 'gzip', compression_opts = 1)

index = 0
for fnam in cxi_in :
    print(fnam)
    with h5py.File(fnam) as f:
        Nevents = f['/entry_1/data_1/data'].shape[0]
        
        with h5py.File(cxi_out, 'r+') as g:
            def write_data(name, obj):
                if type(obj) == h5py.Dataset :
                    if (len(obj.shape) > 0) and (obj.shape[0] == Nevents) :
                        if name not in g :
                            g.create_dataset(name, shape = (frames,) + obj.shape[1:], dtype = obj.dtype, compression = 'gzip', compression_opts = 1, chunks = (1,) + obj.shape[1:])
                        
                        if 'data' not in name :
                            g[name][index: index + len(inds[fnam])] = obj[inds[fnam]]
                        else :
                            #print('writing', name)
                            it =  tqdm(enumerate(inds[fnam]), total = len(inds[fnam]))
                            for i, d in it :
                                it.set_description(f'{name} {index} {i} {d} {index+i}')
                                g[name][index + i] = obj[d]
            
            f.visititems(write_data)
            
    index += len(inds[fnam])


# add soft link 
with h5py.File(cxi_out, 'a') as g:
    g["/entry_1/data_1/data"] = h5py.SoftLink('/entry_1/instrument_1/detector_1/data')

