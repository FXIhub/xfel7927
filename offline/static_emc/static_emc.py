import numpy as np
import h5py
from tqdm import tqdm
import pickle
import os
import sys
import importlib
import runpy

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

from static_emc_init import *

import utils_cl
import utils

# load user defined config.py file
c = utils.load_config(sys.argv[1])

# load recon data 
a = pickle.load(open(c.working_dir + '/recon.pickle', 'rb'))

# i -> pixel index
# d -> frame index
# t -> class index
# l -> background class index

# fluence w[d]
w = a.w

# classes W[t, i]
W = a.W

# background weights for each frame b[d, l]
b = a.b

# background classes B[l, i]
B = a.B

# log likelihood for each frame and class LR[d, c]
LR = a.LR

# probability matrix P[d, c]
P = a.P

# photon counts K[d, i]
K, inds = pickle.load(open(c.working_dir + '/photons.pickle', 'rb'))

if rank == 0 :
    print('classes    :', W.shape[0])
    print('frames     :', len(K))
    print('pixels     :', W.shape[1])
    print('iterations :', a.iterations)


def save_iteration(a):
    # save everything except K and inds
    b = A(a.C, a.L, a.D, a.I, a.mask, a.B, a.pixel_indices, a.file_index, a.frame_index, a.frame_shape, a.frame_slice, a.beta)
    b.C = a.C
    b.L = a.L
    b.D = a.D
    b.I = a.I
    b.frame_shape = a.frame_shape
    b.mask = a.mask
    b.B = a.B
    b.pixel_indices = a.pixel_indices
    b.beta = a.beta
    b.w = a.w
    b.W = a.W
    b.b = a.b
    b.B = a.B
    b.most_likely_classes = a.most_likely_classes
    b.expectation_values = a.expectation_values
    b.LL = a.LL
    b.iterations = a.iterations
    b.P = a.P
    
    fnam = c.working_dir + '/recon_%.4d.pickle'%a.iterations
    print('writing', fnam)
    pickle.dump(b, open(fnam, 'wb'))
    os.system(f"cp {fnam} {c.working_dir}/recon.pickle")
    

for i in range(c.iters):
    beta = c.betas[i]
    update_b = c.update_b[i]
    update_B = c.update_B[i]
    tol_P = c.tol_P
    
    LL, E = utils_cl.calculate_P(K, inds, w, W, b, B, LR, P, beta)
    
    utils_cl.update_w(P, w, W, b, B, K, inds, tol_P = tol_P, min_val = 1e-3, update_b = update_b)
    
    utils_cl.update_W(P, w, W, K, inds, b, B, tol_P = tol_P)
    
    if update_B : utils_cl.update_B(P, w, W, K, inds, b, B, tol_P = tol_P, minval = 1e-10)
    
    # keep track of log-likelihood values
    a.beta = beta
    a.most_likely_classes.append(np.argmax(P, axis=1))
    a.LL.append(LL)
    a.expectation_values.append(E)
    a.iterations += 1
    #utils.plot_iter(a, a.iterations)
    #os.system("pdfunite recon_*.pdf recon.pdf")
    
    # save state
    #if rank == 0 : pickle.dump(a, open('recon.pickle', 'wb'))
    if rank == 0 : save_iteration(a)

