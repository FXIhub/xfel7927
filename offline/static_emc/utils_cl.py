import numpy as np
import pyopencl as cl
import pyopencl.array 
from tqdm import tqdm
import sys
import pathlib

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

gpu_precision = np.float32

# find an opencl device (preferably a GPU) in one of the available platforms
done = False
for p in cl.get_platforms():
    devices = p.get_devices(cl.device_type.GPU)
    if (len(devices) > 0) and ('NVIDIA' in p.name):
        done = True
        break

if not done :
    for p in cl.get_platforms():
        devices = p.get_devices(cl.device_type.GPU)
        if (len(devices) > 0) :
            break
    
if len(devices) == 0 :
    for p in cl.get_platforms():
        devices = p.get_devices()
        if len(devices) > 0:
            break

print('number of devices:', len(devices))
print(rank, 'my device:', devices[rank % len(devices)])
sys.stdout.flush()
context = cl.Context([devices[rank % len(devices)]])
queue   = cl.CommandQueue(context)

code = pathlib.Path(__file__).resolve().parent.joinpath('utils.c')
cl_code = cl.Program(context, open(code, 'r').read()).build()


# make an iterator that splits N into chunks of size n
class chunk_csize:
    def __init__(self, N, n):
        self.chunks = int(np.ceil(N/n))
        self.istart = np.int32(-n)
        self.n      = n
        self.N      = N
        self.counter = 0
    
    def __iter__(self):
        return self

    def __next__(self):
        self.counter += 1
        self.istart  += self.n
        self.istop   = np.int32(min(self.istart + self.n, self.N))
        if self.counter <= self.chunks :
            return (self.istart, self.istop, np.int32(self.istop - self.istart))
        raise StopIteration

    def __len__(self):
        return self.chunks

        

def calculate_LR(K, inds, w, W, b, B, LR, beta, min_val = 1e-10):
    """
    each rank processes all frames for a given set of classes
    LR[d, t] = sum_i K[d, i] log( T[t, d, i] ) - T[t, d, i]
    tranpose for faster summing
    LR[d, t] = sum_i K[i, d] log( T[i, t, d] ) - T[i, t, d]
    """
    D = np.int32(w.shape[0])
    C = np.int32(W.shape[0])
    L = np.int32(b.shape[1])
    I = np.int32(W.shape[1])
    
    beta = np.float32(beta)
    if rank == 0 : print('\nProbability matrix calculation')
    
    # split classes by rank
    my_classes = list(range(rank, C, size))
    classes    = np.int32(len(my_classes))
    
    # parallelise over d chunks on gpu
    frames = np.int32(min(50000, D))
    
    LR_cl = cl.array.zeros(queue, (frames,), dtype = np.float32)
    w_cl  = cl.array.empty(queue, (frames,), dtype = np.float32)
    W_cl  = cl.array.empty(queue, (I,),      dtype = np.float32)
    b_cl  = cl.array.empty(queue, (frames, L),       dtype = np.float32)
    B_cl  = cl.array.empty(queue, (L, I),            dtype = np.float32)
    K_cl  = cl.array.empty(queue, (I, frames),       dtype = np.uint8)
    
    K_dense = np.zeros((I, frames), dtype = np.uint8)
    
    LR_buf = np.empty((frames,), dtype = np.float32)
    
    if rank == 0 :
        disable = False
    else :
        disable = True
    
    # loop over classes
    for c in tqdm(my_classes, desc='processing class', disable = disable):
        cl.enqueue_copy(queue, W_cl.data, W[c])
        cl.enqueue_copy(queue, B_cl.data, B)
            
        for dstart, dstop, dd in tqdm(chunk_csize(D, frames), desc = 'processing frames', leave = False, disable = disable):
            # load transposed photons to gpu 
            K_dense.fill(0)
            for i, d in enumerate(range(dstart, dstop)):
                K_dense[inds[d], i] = K[d]
            cl.enqueue_copy(queue, K_cl.data, K_dense)
            
            # load class etc. to gpu
            cl.enqueue_copy(queue, w_cl.data, w[dstart:dstop])
            cl.enqueue_copy(queue, b_cl.data, b[dstart:dstop])
            
            cl_code.calculate_LR_T_dt(queue, (dd,), None, 
                    LR_cl.data,  K_cl.data, w_cl.data, W_cl.data,
                    b_cl.data, B_cl.data, beta, L, I, frames)
            
            sys.stdout.flush()
            cl.enqueue_copy(queue, LR_buf, LR_cl.data)
            LR[dstart: dstop, c] = LR_buf[:dd]
    
    # all gather
    #print('gathering LR')
    #sys.stdout.flush()
    for r in range(size):
        #print(size, rank, r)
        #sys.stdout.flush()
        rank_classes = list(range(r, C, size))
        LR[:, rank_classes] = comm.bcast(LR[:, rank_classes], root=r)
    
def calculate_expectation(K, inds, w, W, b, B, LR, P, beta):
    calculate_LR(K, inds, w, W, b, B, LR, beta)
    
    # E = sum_dt P[d, t] LR[d, t]
    expectation = np.sum(P * LR)
    return expectation

def calculate_P(K, inds, w, W, b, B, LR, P, beta):
    LR.fill(0)
    calculate_LR(K, inds, w, W, b, B, LR, beta)
    
    #print('Allreducing probablility matrix') 
    #sys.stdout.flush()
    #x = np.zeros_like(LR)
    #comm.Allreduce(LR, x, op = MPI.SUM)
    #LR[:] = x

    
    # calculate log-likelihood before normalisation
    LL = np.sum(LR)
    
    # E = sum_dt P[d, t] LR[d, t]
    expectation = np.sum(P * LR) / beta
    
    normalise_LR(LR, P)
    return LL, expectation
    
def normalise_LR(LR, P):
    """
    P[d,t] = exp(LR[d,t]) / sum_t exp(LR[d,t])
    """
    # normalise to produce probability matrix
    m = np.max(LR, axis=1)
    
    LR[:] = np.exp(LR - m[:, None])
    P[:]  = LR / np.sum(LR, axis=-1)[:, None]
    
def update_W(P, w, W, K, inds, b, B, tol_P = 1e-2, minval = 1e-10, update_B = True):
    """
    use dense K so we can parallelise
    but still mask based on P[t, d]
    
    grad[t, i] = sum_d P[d, t] w[d] ( K[i, d] / T[t, d, i] - 1)
               = sum_d P[d, t] w[d] K[i, d] / (w[d] W[t, i] + B[d, i]) - sum_d P[d, t] w[d] 
               = sum_d P[d, t] K[i, d] / (W[t, i] + B[d, i] / w[d]) - g0[t]
    
    - fast (sum over slow axis)
    g0[t] = sum_d P[d, t] w[d]
    
    - fast (d mask + sum over slow axis)
    - but we will only be able to store K[d, i] when P is sparse
    - then we will have to chunk over d 
    - one worker per t i index
    xmax[t, i] = sum_d P[d, t] K[d, i] / g0[t]
    
    - fast
    f[t, i] =  sum_d P[d, t] K[d, i] / (W[t, i] + B[d, i] / w[d]) 
    g[t, i] = -sum_d P[d, t] K[d, i] / (W[t, i] + B[d, i] / w[d])^2 
    
    - do not want to store on gpu
    B[d, i] = b[d, l] B[l, i]
    """
    C = np.int32(W.shape[0])
    I = np.int32(W.shape[1])
    D = np.int32(w.shape[0])
    L = np.int32(B.shape[0])
        
    minval = np.float32(minval)
     
    if rank == 0 : print('\nW-update')
    
    my_classes = list(range(rank, C, size))

    # chunk over pixels
    pixels  = np.int32(2048)
    
    K_dense = np.zeros((D, pixels), dtype = np.uint8)
    W_buf   = np.empty((pixels,), dtype = np.float32)
    k       = np.empty((I,), dtype = np.uint8)
    
    K_cl    = cl.array.empty(queue, (D * pixels,), dtype = np.uint8)
    P_cl    = cl.array.empty(queue, (D,), dtype = np.float32)
    w_cl    = cl.array.empty(queue, (D,), dtype = np.float32)
    b_cl    = cl.array.empty(queue, (L * D,), dtype = np.float32)
    B_cl    = cl.array.empty(queue, (L * pixels,), dtype = np.float32)
    W_cl    = cl.array.empty(queue, (pixels,), dtype = np.float32)
    
    if rank == 0 :
        disable = False
    else :
        disable = True
    
    for c in tqdm(my_classes, desc = 'processing class', disable = disable):
        # mask low P frames
        p      = P[:, c]
        frames = np.where(p > (p.max() * tol_P))[0]
        Dc     = np.int32(len(frames))
        Pd     = np.ascontiguousarray(P[frames, c])
        wd     = np.ascontiguousarray(w[frames])
        bdl    = np.ascontiguousarray(b[frames])

        # calculate g0
        gW = np.float32(np.sum(Pd * wd))
        
        # load P
        cl.enqueue_copy(queue, P_cl.data, Pd[:Dc])
        
        # load w
        cl.enqueue_copy(queue, w_cl.data, wd)

        # load b
        cl.enqueue_copy(queue, b_cl.data, bdl)
            
        for istart, istop, di in tqdm(chunk_csize(I, pixels), desc = 'processing pixels', leave = False, disable = disable):
            # load K
            for i, d in enumerate(frames) :
                k.fill(0)
                k[inds[d]] = K[d]
                K_dense[i][:di] = k[istart:istop]
            
            cl.enqueue_copy(queue, K_cl.data, K_dense[:Dc])
            
            # load background to gpu 
            cl.enqueue_copy(queue, B_cl.data, np.ascontiguousarray(B[:, istart:istop]))
            
            # load W
            cl.enqueue_copy(queue, W_cl.data, np.ascontiguousarray(W[c, istart:istop]))
            
            cl_code.update_W(queue, (di,), None, 
                        P_cl.data, K_cl.data, b_cl.data, B_cl.data, 
                        w_cl.data, W_cl.data, gW, minval, pixels, L, Dc)
            
            cl.enqueue_copy(queue, W_buf, W_cl.data)
            W[c, istart:istop] = W_buf[:di]
    
    # all gather
    for r in range(size):
        rank_classes = list(range(r, C, size))
        W[rank_classes] = comm.bcast(np.clip(W[rank_classes], minval, None), root=r)

def update_B(P, w, W, K, inds, b, B, tol_P = 1e-2, minval = 1e-10, update_B = True):
    """
    grad[l, i] = sum_dt P[d, t] b[d, l] (K[d, i] / T[l, d, t, i] - 1)
               = sum_dt P[d, t] b[d, l] K[d, i] / T[l, d, t, i] - sum_dt P[d, t] b[d, l]
               = sum_dt P[d, t] b[d, l] K[d, i] / (w[d] W[t, i] + sum_l b[d, l] B[l, i]) - sum_d b[d, l]
               = sum_dt P[d, t] b[d, l] K[d, i] / (b[d] B[i] + w[d] W[t, i] + sum_l'neql b[d, l'] B[l', i]) - g0[l]
               = sum_dt P[d, t] K[d, i] / ( B[i] + (w[d] W[t, i] + sum_l'neql b[d, l'] B[l', i]) / b[d]) - g0[l]
    
    g0[l] = sum_t (sum_d P[t, d] b[d, l])
    """
    C = np.int32(W.shape[0])
    I = np.int32(W.shape[1])
    D = np.int32(w.shape[0])
    L = np.int32(B.shape[0])
    
    if rank == 0 : print('\nB-update')

    my_frames = list(range(rank, D, size))
    
    if rank == 0 :
        disable = False
    else :
        disable = True

    # calculate xmax[l, i] = sum_dt P[d, t] K[d, i] / sum_dt P[d, t] b[d, l]
    #                      = sum_d K[d, i] / sum_d b[d, l]
    
    Ksums = np.zeros((I,), dtype = int)
    for d in range(D):
        Ksums[inds[d]] += K[d]
    
    bsums = np.sum(b, axis = 0)
    
    g = np.empty((I,), dtype = float)
    f = np.empty((I,), dtype = float)
    g_buf = np.empty((I,), dtype = float)
    f_buf = np.empty((I,), dtype = float)
    for l in tqdm(range(L), desc = 'processing class', disable = disable):
        # xmax[i]
        xmax = Ksums / bsums[l]
        
        x  = B[l].copy()
        c  = bsums[l]
        
        ls = [ll for ll in range(L) if ll != l]
        b2 = b[:, ls] 
        B2 = B[ls] 

        for iter in range(3):
            g_buf.fill(0)
            f_buf.fill(0)
            # sum over d to calculate g and f
            for d in tqdm(my_frames, desc = 'processing frames', disable = disable, leave = False):
                p  = P[d] 
                ts = np.where(p > (p.max() * tol_P))[0]
                t, i = np.ix_(ts, inds[d])
                
                # could gpu accellerate this, but is it worth it?
                # could be, since there may be a lot of t's
                # but I don't plan to update B for large datasets
                if L > 1 :
                    T  = x[i] + (w[d] * W[t, i] + np.dot(b2[d], B2[:, inds[d]])) / b[d, l]
                else :
                    T  = x[i] + (w[d] * W[t, i]) / b[d, l]
                
                PK = P[d, t] * K[d]
                
                f_buf[i]  += np.sum(PK / T,    axis=0)
                g_buf[i]  -= np.sum(PK / T**2, axis=0)

            f.fill(0)
            g.fill(0)
            comm.Allreduce(f_buf, f, op = MPI.SUM)
            comm.Allreduce(g_buf, g, op = MPI.SUM)
            
            u  = - f * f / g 
            v  = - f / g - x 
            xp = u / c - v 
             
            x = np.clip(xp, minval, xmax) 
        
        B[l] = x
    
            
def update_w(P, w, W, b, B, K, inds, tol_P = 1e-3, tol = 1e-5, min_val = 1e-10, max_iters=1000, update_b = True):
    """
    gw[d] = sum_t P[d, t] sum_i W[t, i] (K[d, i] / T[i, t, d] - 1)
          = sum_t P[d, t] sum_i W[t, i] K[d, i] / T[i, t, d] - sum_t P[d, t] sum_i W[t, i]
          = sum_t P[d, t] sum_i W[t, i] K[d, i] / (w[d] W[t, i] + sum_l b[d, l] B[l, i]) - sum_t P[d, t] sum_i W[t, i]
          = sum_t P[d, t] sum_i K[d, i] / (w[d] + sum_l b[d, l] B[l, i] / W[t, i]) - g0[d]
    
    
    gb[d, l] = sum_t P[d, t] sum_i B[l, i] (K[d, i] / T[i, t, d, l] - 1)
             = sum_t P[d, t] sum_i B[l, i] K[d, i] / (w[d] W[t, i] + sum_l' b[d, l'] B[l', i]) - (sum_t P[d, t]) (sum_i B[l, i])
             = sum_t P[d, t] sum_i K[d, i] / (b[d, l] + (w[d] W[t, i] + sum_l'!=l b[d, l'] B[l', i])/B[l, i]) - sum_i B[l, i]
    """
    # check if things will fit in memory
    C = np.int32(W.shape[0])
    I = np.int32(W.shape[1])
    D = np.int32(w.shape[0])
    L = np.int32(B.shape[0])
    
    # each rank processes a sub-set of frames
    my_frames = list(range(rank, D, size))
    
    if rank == 0 :
        print('\nupdating fluence estimates')
        disable = False
    else :
        disable = True
    
    Wsums = np.sum(W, axis=-1)
    g0    = np.dot(P[my_frames], Wsums)
    Ksums = np.array([np.sum(K[d]) for d in my_frames])
    xmax  = Ksums / g0
    
    if update_b :
        gb0    = np.sum(B, axis = -1)
        xmax_b = Ksums[:, None] / gb0[None, :]
    
    for j, d in tqdm(enumerate(my_frames), total = len(my_frames), desc = 'processing frame', disable = disable):
        p = P[d] 
        classes    = np.where(p > (p.max() * tol_P))[0]
        t, i       = np.ix_(classes, inds[d])
        
        Wti  = np.ascontiguousarray(W[t, i])
        Bi   = np.dot(b[d], B[:, inds[d]]) / Wti
        Pt   = np.ascontiguousarray(P[d][t])
        PK   = np.ascontiguousarray(Pt * K[d])
        x    = np.clip(w[d], min_val, xmax[j])
        c    = g0[j]

        # w update
        for iter in range(5):
            T = x + Bi
            f =  np.sum(PK / T)
            g = -np.sum(PK / T**2)
             
            u  = - f * f / g 
            v  = - f / g - x 
            xp = u / c - v 
             
            x = np.clip(xp, min_val, xmax[j]) 
        
        w[d] = x;
        
        # gb[d, l] = sum_t P[d, t] sum_i K[d, i] / (b[d, l] + (w[d] W[t, i] + sum_l'!=l b[d, l'] B[l', i])/B[l, i]) - sum_i B[l, i]
        # b update
        if update_b :
            for l in range(L):
                ls = [ll for ll in range(L) if ll != l]
                c = gb0[l]
                x = np.clip(b[d, l], min_val, xmax_b[j, l])
                if L > 1 :
                    B2i = (w[d] * Wti + np.dot(b[d, ls], B[np.ix_(ls, inds[d])])) / B[l][inds[d]]
                else :
                    B2i = w[d] * Wti / B[l, inds[d]]
                  
                for iter in range(5):
                    T = x + B2i 
                    f =  np.sum(PK / T)
                    g = -np.sum(PK / T**2)
                     
                    u  = - f * f / g 
                    v  = - f / g - x 
                    xp = u / c - v 
                     
                    x = np.clip(xp, min_val, xmax_b[j, l]) 
                
                b[d, l] = x;
        
    # all gather
    for r in range(size):
        rank_frames = list(range(r, D, size))
        w[rank_frames] = comm.bcast(w[rank_frames], root=r)
    
        if update_b : b[rank_frames, :] = comm.bcast(b[rank_frames, :], root=r)
    


