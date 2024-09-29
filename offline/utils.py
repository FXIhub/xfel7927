import numpy as np
import h5py

    
def update_h5(f, key, value, compression = False):
    if key in f and f[key].shape != value.shape:
        del f[key]
        
    if key in f and f[key].shape == value.shape:
        f[key][...] = value
    
    if key not in f:
        if compression :
            f.create_dataset(key, data = value, compression='gzip', compression_opts=1, shuffle=True, chunks=True)
        else :
            f.create_dataset(key, data = value)

def put_a_in_b_by_train(a, a_tr, a_ce, b, b_tr, b_ce):
    # align by train then pulse
    trains = np.unique(a_tr)
    a_dict = {t:{} for t in trains} 
    for v, train, cell in zip(a, a_tr, a_ce):
        a_dict[train][cell] = v
    
    for i, train, cell in zip(range(b.shape[0]), b_tr, b_ce):
        if train in a_dict and cell in a_dict[train]:
            b[i] = a_dict[train][cell]
    return b
        
        

prefix = {'y': 1e-24,  # yocto
          'z': 1e-21,  # zepto
          'a': 1e-18,  # atto
          'f': 1e-15,  # femto
          'p': 1e-12,  # pico
          'n': 1e-9,   # nano
          'u': 1e-6,   # micro
          'm': 1e-3,   # mili
          'c': 1e-2,   # centi
          'd': 1e-1,   # deci
          'k': 1e3,    # kilo
          'M': 1e6,    # mega
          'G': 1e9,    # giga
          'T': 1e12,   # tera
          'P': 1e15,   # peta
          'E': 1e18,   # exa
          'Z': 1e21,   # zetta
          'Y': 1e24,   # yotta
    }
