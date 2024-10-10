import numpy as np

PREFIX = '/home/andyofmelbourne/Documents/2023/P3004-take-2/gold/'
data = []
for i in range(87, 96):
    data.append(PREFIX + f'hits_r00{i}.cxi')

classes            = 16
rotations          = 64
sampling           = 4
model_length       = 64
background_classes = 1
max_frames         = 8192
#frame_shape        = (16, 128, 512)
#frame_slice        = np.s_[:, :, :]
#imshow             = lambda x: x[frame_slice]
#pixels             = imshow(np.arange(1024**2).reshape(frame_shape))
polarisation_axis  = 'x'
filter_by          = '/static_emc/good_hit'
filter_value       = True
dtype              = np.float32


tol_P = 1e-2

iters = 3
update_b = np.ones((iters,), dtype=bool)
update_b[0] = False
update_B = np.zeros((iters,), dtype=bool)
betas = iters*[0.001]
no_back = iters*[False]
#beta_start = 0.001
#beta_stop  = 0.01
#betas = (beta_stop / beta_start)**(np.arange(iters)/(iters-1)) * beta_start
