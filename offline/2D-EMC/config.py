import numpy as np

data = ['/home/andyofmelbourne/Documents/2024/p7927/scratch/saved_hits/Ery_size_filtered.cxi']

classes            = 16
rotations          = 64
sampling           = 4
model_length       = 64
background_classes = 1
max_frames         = 2048
#frame_shape        = (16, 128, 512)
#frame_slice        = np.s_[:, :, :]
#imshow             = lambda x: x[frame_slice]
#pixels             = imshow(np.arange(1024**2).reshape(frame_shape))
polarisation_axis  = 'x'
filter_by          = None # '/entry_1/instrument_1/detector_1/hit_sigma' #'/static_emc/good_hit'
filter_value       = 100 # True
dtype              = np.float32


tol_P = 1e-2

iters = 3
update_b = np.ones((iters,), dtype=bool)
update_b[:] = False
update_B = np.zeros((iters,), dtype=bool)
betas = iters*[0.2]
no_back = iters*[True]
#beta_start = 0.001
#beta_stop  = 0.01
#betas = (beta_stop / beta_start)**(np.arange(iters)/(iters-1)) * beta_start
