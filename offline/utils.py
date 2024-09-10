import numpy as np
import h5py

    
def update_h5(f, key, value):
    if key in f and f[key].shape != value.shape:
        del f[key]
        
    if key in f and f[key].shape == value.shape:
        f[key][...] = value
    
    if key not in f:
        f[key] = value
