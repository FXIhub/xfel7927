import numpy as np
import h5py
import pickle


I = pickle.load(open('merged_intensity.pickle', 'rb'))['I']


with h5py.File('merged_intensity.h5', 'w') as f:
    f.create_dataset('I', data = I, dtype = np.float32, compression = 'gzip')
