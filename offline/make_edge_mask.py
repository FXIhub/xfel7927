import numpy as np
import h5py

def make_edge_mask(shape = (512, 128)):
    edges_mask = np.ones(shape, dtype = bool)
    edges_mask[::64] = False
    edges_mask[63::64] = False
    edges_mask[:, 0] = False
    edges_mask[:, -1] = False
    return edges_mask

shape = (16, 512, 128)

mask = np.ones(shape, dtype=bool)
edge_mask_panel = make_edge_mask()

for i in range(mask.shape[0]):
    mask[i] = edge_mask_panel

with h5py.File('panel_edge_mask.h5', 'w') as f:
    f['data'] = mask

    # optional
    #f['pixel_size'] = 200e-6
