import numpy as np
import extra_geom
import h5py

# to be filled in by script
data = []
working_dir = ''
static_emc_mask_fnam = ''
geom_fnam = ''

geom = extra_geom.AGIPD_1MGeometry.from_crystfel_geom(geom_fnam)

with h5py.File(static_emc_mask_fnam) as f:
    static_emc_mask = f['entry_1/good_pixels'][()]

def crop_to_nonzero_ND(ar):
    axes = list(range(ar.ndim))
    mins  = []
    maxes = []
    for d in axes:
        a = list(axes)
        a.pop(d)
        s = np.nansum(ar, axis = tuple(a))
        mins.append(np.where(s)[0][0])
        maxes.append(ar.shape[d] - np.where(s[::-1])[0][0])

    sl = tuple(slice(i, j) for i, j in zip(mins,maxes))
    return ar[sl].copy(), sl


classes            = 100
background_classes = 1
max_frames         = 100000
frame_shape        = (16, 512, 128)
#frame_slice        = np.s_[4, -128:, -128:]
#imshow             = lambda x: x[frame_slice]
frame_slice        = static_emc_mask
imshow             = lambda x: crop_to_nonzero_ND(geom.position_modules(x)[0])[0]
# just use the first part of the first panel (low q)
#pixels             = imshow(np.arange(1024**2).reshape(frame_shape))
filter_by          = None
filter_value       = True

tol_P = 1e-2

iters = 20
update_b     = np.ones((iters,), dtype=bool)
update_b[0]  = False
update_B     = np.zeros((iters,), dtype=bool)
beta_start   = 0.001
beta_stop    = 0.1
betas = (beta_stop / beta_start)**(np.arange(iters)/(iters-1)) * beta_start

