import numpy as np
import scipy.constants as sc
#wav = sc.h * sc.c / photon_energy
from scipy.interpolate import RegularGridInterpolator
from tqdm import tqdm
import pathlib
import runpy
import h5py


try :
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0 :
        silent = False
    else :
        silent = True
except ImportError :
    silent = False

def load_config(path):
    p = pathlib.Path(path)
    
    # returns a dict
    config = runpy.run_path(str(p.absolute()))
    
    return config


def plot_iter(recon_dir):
    """
    | P-matrix |
    | lines P  |
    | w, b     |
    | W most   |
    | W middle |
    | W least  |
    | favour   |
    | B        | 
    | LL       |
    """
    import matplotlib.pyplot as plt

    recon_file = recon_dir + '/recon.h5'
    data_file = recon_dir + '/data.cxi'
    
    layout = """
        PPP012
        PPP345
        www678
        LLLfff
    """
    #fig = plt.figure(constrained_layout=True)
    fig = plt.figure(tight_layout=True)
    fig.set_size_inches(20, 15)
    ax_dict = fig.subplot_mosaic(layout)

    with h5py.File(recon_file) as f:
        most_likely_classes = f['/iterations/most_likely_class'][()]
        w                   = f['fluence'][()]
        b                   = f['background_weights'][()]
        P                   = f['probability_matrix'][()]
        I                   = f['models'][()]
        B                   = f['background'][()]
        expectation_values  = f['/iterations/expectation_value'][()]
        iteration           = f['/iterations/iters'][()] - 1

    with h5py.File(data_file) as f:
        pixel_indices = f['/entry_1/instrument_1/detector_1/pixel_indices'][()]
        frame_shape = f['/entry_1/instrument_1/detector_1/frame_shape'][()]
    
    D = most_likely_classes.shape[1]

    c = load_config(recon_dir + '/config.py')
    
    # P-matrix
    ax = ax_dict["P"]
    # assign colour (value) to each class in subset
    frames  = np.random.randint(0, D, size = min(256, D))
    array   = np.zeros((iteration+1, frames.shape[0]), dtype=float)
    array[-1, :] = np.linspace(0, 1, frames.shape[0])
    for i in range(iteration-1, -1, -1):
        a1 = most_likely_classes[i][frames]
        a2 = most_likely_classes[i+1][frames]
        for j in range(array.shape[1]):
            if (a1[j] - a2[j]) != 0 :
                array[i, j] = np.random.random()
            else :
                array[i, j] = array[i+1, j]
    im = ax.imshow(array.T, aspect = 'auto', origin='lower', interpolation='nearest')#, cmap = 'gist_ncar_r')
    ax.set_ylabel('frame')
    ax.set_xlabel('iteration')
    ax.set_title('most likely class')
    ax.set_xticks(range(iteration+1))
    
    # w
    ax = ax_dict["w"]
    ax.plot(range(w.shape[0]), w, linewidth=0.8, alpha=1.)
    ax.set_xlabel('frame number (d)')
    ax.set_ylabel('relative fluence w[d]')
    ax.spines[['top']].set_visible(False)
    
    # b
    ax2 = ax.twinx()
    for l in range(b.shape[1]):
        ax2.plot(range(b.shape[0]), b[:, l], linewidth=0.8, alpha=.7, c = 'k', label=f'b class {l}')
    ax2.set_ylabel('background weight')
    ax2.spines[['top']].set_visible(False)

    # W 
    # favour = sum_d P[d, t, r]
    C = I.shape[0]
    favour = np.sum(P, axis=(0, 2))
    ts     = np.argsort(favour)
    most   = ts[-3:][::-1]
    middle = ts[C//2-1:C//2 + 2][::-1]
    least  = ts[:3][::-1]
    
    fav     = [favour[most], favour[middle], favour[least]]
    labels  = [most, middle, least]
    classes = [I[most], I[middle], I[least]]
    
    # show 3 most favoured classes and 3 least favoured class
    for i in range(3):
        for j in range(3):
            k = str(3 * i + j)
            ax = ax_dict[k]
            ax.imshow(classes[i][j]**0.2, origin='lower')
            ax.axis('off')
            ax.set_title(f'class {labels[i][j]} no. of frames {round(fav[i][j])}', fontsize=5)

    # favour bar plot
    ax = ax_dict["f"]
    ax.bar(np.arange(C), favour[ts[::-1]],  width = 1, align='edge', color='lightcoral', edgecolor='k', alpha=0.8, linewidth=1)
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_xlim([0, C])
    ax.set_title("class favour")
    ax.set_ylabel("occupancy")
    ax.set_xlabel("classes (sorted)")
    
    # expectation value plot
    ax = ax_dict["L"]
    ax.bar(np.arange(len(expectation_values)), expectation_values,  width = 1, align='edge', color='lightcoral', edgecolor='k', alpha=0.8, linewidth=1)
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_xlim([0, max(100, len(expectation_values))])
    #ax.set_yscale('log')
    ax.set_title("expectation values")
    ax.set_xlabel("iterations")
    
    print('saving', recon_dir + f'/recon_{str(iteration).zfill(3)}.pdf')
    plt.savefig(recon_dir + f'/recon_{str(iteration).zfill(3)}.pdf')

    plt.close(fig)

def save(outdir, w, b, P, logR, I, beta, expectation_value, log_likihood, iteration):
    with h5py.File(outdir + '/recon.h5', 'a') as f:
        f['models'][...] = I
        f['fluence'][...] = w
        f['background_weights'][...] = b
        f['probability_matrix'][...] = P
        f['logR'][...] = logR
        
        # resize datasets if required
        if f['iterations/expectation_value'].shape[0] <= iteration :
            f['iterations/expectation_value'].resize(iteration + 1, axis=0)
            f['iterations/log_likihood'].resize(iteration + 1, axis=0)
            f['iterations/most_likely_class'].resize(iteration + 1, axis=0)
            f['iterations/beta'].resize(iteration + 1, axis=0)
        
        f['iterations/expectation_value'][iteration] = expectation_value
        f['iterations/log_likihood'][iteration] = log_likihood
        f['iterations/beta'][iteration] = beta
        mc = np.argmax(np.max(P, axis=2), axis=1)
        f['iterations/most_likely_class'][iteration] = mc
        f['iterations/iters'][...] = iteration + 1
        
        
            

def write_data_file(output, fnams, max_frames, frame_slice, filter_by, filter_value, polarisation_axis = 'x'):
    # write meta data
    with h5py.File(fnams[0], 'r') as f:
        data_dtype  = f['entry_1/data_1/data'].dtype
        frame_shape = f['entry_1/data_1/data'].shape[1:]
        mask_temp   = f['/entry_1/instrument_1/detector_1/mask'][()]
        mask = np.zeros_like(mask_temp)
        mask[frame_slice] = mask_temp[frame_slice]
        
        datasets = ['/entry_1/instrument_1/detector_1/mask',
                    '/entry_1/instrument_1/detector_1/xyz_map',
                    '/entry_1/instrument_1/detector_1/x_pixel_size',
                    '/entry_1/instrument_1/detector_1/y_pixel_size',
                    '/entry_1/instrument_1/detector_1/pixel_area',
                    '/entry_1/instrument_1/detector_1/distance',
                    '/entry_1/instrument_1/detector_1/background',
                    '/entry_1/instrument_1/source_1/energy']
        
        file_metadata = {}
        for d in datasets:
            v = f[d][()]
            
            # assume any dataset with frame shape as the last dimensions
            # need to be masked (including mask)
            if v.shape[-len(frame_shape):] == frame_shape :
                v = v[..., mask]
            
            file_metadata[d] = v.copy()
        
        file_metadata['/entry_1/instrument_1/detector_1/pixel_indices'] = np.where(mask.ravel())[0]
        
        with h5py.File(output, 'w') as g:
            for d, v in file_metadata.items():
                print('writing dataset', d,'from first input file')
                if d in g :
                    g[d][...] = v
                else :
                    g[d] = v
    
    # get location of good frames in datasets
    frames = {}
    number_of_frames = 0 
    for fnam in fnams :
        with h5py.File(fnam, 'r') as f:
            frames[fnam] = np.where(f[filter_by][()] == filter_value)[0]
            number_of_frames += len(frames[fnam])
            
            if number_of_frames >= max_frames :
                frames[fnam] = frames[fnam][:max_frames - number_of_frames]
                number_of_frames = max_frames
                break
    
    print('found', number_of_frames, 'to write to output file')
    file_metadata['frames'] = number_of_frames
    
    out_shape = (np.sum(mask),)
    
    frame_index = 0
    with h5py.File(output, 'a') as g:
        data = g.create_dataset('entry_1/data_1/data', 
                                (number_of_frames,) + out_shape,  
                                dtype = data_dtype, 
                                compression='gzip',
                                compression_opts = 1)
        
        for fnam in frames :
            
            with h5py.File(fnam, 'r') as f:
                data_file = f['entry_1/data_1/data']
                
                for d in tqdm(frames[fnam], desc = fnam):
                    frame             = data_file[d][mask] 
                    data[frame_index] = frame
                    frame_index      += 1
    return file_metadata


def solid_angle_polarisation_factor(xyz, pixel_area, polarisation_axis = 'x'):
    
    C      = np.zeros(xyz.shape[1:], dtype = float)
    radius = np.sum(xyz**2, axis = 0)**0.5

    # polarisation correction
    if polarisation_axis == 'x' :
        C[:] = 1 - (xyz[0] / radius)**2 
    elif polarisation_axis == 'y' :
        C[:] = 1 - (xyz[1] / radius)**2 
    else :
        raise ValueError("polarisation axis must be one of 'x' or 'y'")

    # solid angle correction
    C[:] *= pixel_area * xyz[2] / radius**3 

    # rescale
    C /= C.max()
    return C


def calculate_rotation_matrices(M):
    # theta[r] = 2 pi r / M_in_plane
    # R[r]     = [cos -sin]
    #            |sin  cos|
    t = 2 * np.pi * np.arange(M) / M
    R = np.empty((M, 2, 2), dtype = np.float32)
    R[:, 0, 0] =  np.cos(t)
    R[:, 0, 1] = -np.sin(t)
    R[:, 1, 0] =  np.sin(t)
    R[:, 1, 1] =  np.cos(t)
    return R

def expand(classes, points_I, I, W, xyz, R, C, minval = 1e-8):
    _, rotations, pixels = W.shape
    
    interp = RegularGridInterpolator(points_I, I[0], fill_value = 0.)
    points = np.empty((pixels, 2))
    for c in tqdm(classes, desc = 'Expand', disable = silent):
        interp.values[:] = I[c]
        for r in range(rotations):
            points[:, 0] = R[r, 0, 0] * xyz[0] + R[r, 0, 1] * xyz[1]
            points[:, 1] = R[r, 1, 0] * xyz[0] + R[r, 1, 1] * xyz[1]
            W[c, r] = np.clip(C * interp(points), minval, None)

def compress(classes, P, K, W, R, xyz, i0, dx, I):
    _, rotations, pixels = W.shape
    
    I.fill(0)
    overlap1 = np.zeros(I.shape, dtype=int)
    overlap2 = np.zeros(I.shape, dtype=float)
    #weights = np.tensordot(np.sum(K, axis=1), P, axes=1)
    weights = np.sum(P, axis=0)
    print(weights)
    
    points = np.empty((2, pixels))
    mi     = np.empty((pixels,), dtype=int)
    mj     = np.empty((pixels,), dtype=int)
    for c in tqdm(classes, desc = 'Compress', disable = silent):
        for r in range(rotations):
            points[0, :] = R[r, 0, 0] * xyz[0] + R[r, 0, 1] * xyz[1]
            points[1, :] = R[r, 1, 0] * xyz[0] + R[r, 1, 1] * xyz[1]
            mi[:]   = np.round(i0 + points[0]/dx)
            mj[:]   = np.round(i0 + points[1]/dx)
            
            np.add.at(I[c], (mi, mj), W[c, r] * weights[c, r])
            np.add.at(overlap1[c], (mi, mj), 1)
            np.add.at(overlap2[c], (mi, mj), weights[c, r])

    overlap2[overlap1 == 0] = 1
    I /= overlap2

def compress_P_weight(classes, P, W, R, xyz, i0, dx, I, tol_P = 1e-3):
    _, rotations, pixels = W.shape
    
    I.fill(0)
    overlap1 = np.zeros(I.shape, dtype=int)
    overlap2 = np.zeros(I.shape, dtype=float)
    #weights = np.tensordot(np.sum(K, axis=1), P, axes=1)
    weights = np.sum(P, axis=0)
    
    points = np.empty((2, pixels))
    mi     = np.empty((pixels,), dtype=int)
    mj     = np.empty((pixels,), dtype=int)
    for c in tqdm(classes, desc = 'Compress', disable = silent):
        for r in np.where(weights[c] > (tol_P * weights[c].max()))[0]:
            points[0, :] = R[r, 0, 0] * xyz[0] + R[r, 0, 1] * xyz[1]
            points[1, :] = R[r, 1, 0] * xyz[0] + R[r, 1, 1] * xyz[1]
            mi[:]   = np.round(i0 + points[0]/dx)
            mj[:]   = np.round(i0 + points[1]/dx)
            
            np.add.at(I[c], (mi, mj), W[c, r] * weights[c, r])
            np.add.at(overlap1[c], (mi, mj), 1)
            np.add.at(overlap2[c], (mi, mj), weights[c, r])
    
    overlap2[overlap1 == 0] = 1
    I /= overlap2

def compress_P_thresh(classes, P, W, R, xyz, i0, dx, I, tol_P = 1e-2):
    _, rotations, pixels = W.shape
    
    I.fill(0)
    overlap1 = np.zeros(I.shape, dtype=int)
    weights  = np.sum(P, axis=0)
    
    points = np.empty((2, pixels))
    mi     = np.empty((pixels,), dtype=int)
    mj     = np.empty((pixels,), dtype=int)
    for c in tqdm(classes, desc = 'Compress', disable = silent):
        for r in np.where(weights[c] > (tol_P * weights[c].max()))[0]:
            points[0, :] = R[r, 0, 0] * xyz[0] + R[r, 0, 1] * xyz[1]
            points[1, :] = R[r, 1, 0] * xyz[0] + R[r, 1, 1] * xyz[1]
            mi[:]   = np.round(i0 + points[0]/dx)
            mj[:]   = np.round(i0 + points[1]/dx)
            
            np.add.at(I[c], (mi, mj), W[c, r] * weights[c, r])
            np.add.at(overlap1[c], (mi, mj), 1)

    overlap1[overlap1 == 0] = 1
    I /= overlap1


def calculate_probability_matrix(frames, w, W, I, b, B, K, logR, P, beta):
    # probability matrix
    for i, d in tqdm(enumerate(frames), total = len(frames), desc = 'calculating probability matrix', disable = silent):
        T = w[d] * W + np.dot(b[d], B)
        logR[d] = beta * np.sum(K[d] * np.log(T) - T, axis=-1)
        
        m = np.max(logR[d])
        P[d]     = logR[d] - m
        P[d]     = np.exp(P[d])
        P[d]    /= np.sum(P[d])

def calculate_P_stuff(P, logR, beta):
    expectation_value = np.sum(P * logR) / beta
    log_likihood      = np.sum(logR) / beta
    print('expectation value: {:.2e}'.format(expectation_value))
    print('log likelihood   : {:.2e}'.format(log_likihood))
    return expectation_value, log_likihood

def update_W(my_classes, w, W, b, B, P, K, minval = 1e-15, iters = 4):
    _, rotations, pixels = W.shape
    frames = len(K)
    
    classes = len(my_classes)
    
    # xmax[t, r, i] = sum_d P[d, t, r] K[d, i] / \sum_d w[d] P[d, t, r] 
    c    = np.tensordot(w, P[:, my_classes], axes = ((0,), (0,)))
    c[c<minval] = minval
    
    xmax = np.tensordot(P[:, my_classes], K, axes = ((0,), (0,))) / c[..., None]
    xmax[xmax<minval] = minval
    
    step_min = -xmax/2
    
    #- gW[t, r, i] = sum_d w[d] P[d, t, r] (K[d, i] / T[d, t, r, i] - 1)
    #    = sum_d P[d, t, r] K[d, i] w[d] / T[d, t, r, i]  - \sum_d w[d] P[d, t, r] 
    f  = np.empty((classes, rotations, pixels))
    g  = np.empty((classes, rotations, pixels))
    T  = np.empty((classes, rotations, pixels))
    PK = np.empty((classes, rotations, pixels))
    step = np.empty((classes, rotations, pixels))
    for i in range(iters) :
        f.fill(0)
        g.fill(0)
        for d in tqdm(range(frames), desc = f'updating classes, iteration {i}', disable = silent):
            T[:]  = W[my_classes] + np.dot(b[d], B) / w[d]
            PK[:] = P[d, my_classes, :, None] * K[d]
            f    += PK / T
            g    -= PK / T**2
        
        #print('T^2 min:', (T**2).min())
        #print('T min:', T.min())
        #print('g min max:', g.min(), g.max())
        #print('f min max:', f.min(), f.max())
        step[:] = f / g * (1 - f / c[..., None])
        #print('step min max:', step.min(), step.max())
        #print('c min:', c.min())
        g[g > -minval] = -minval
        step[:] = f / g * (1 - f / c[..., None])
        np.clip(step, step_min, None, step)
        W[my_classes]   += step
        m = W[my_classes] > minval
        W[my_classes] = np.clip(W[my_classes], minval, xmax)    
        
        if not silent : 
            print(i, np.mean((f - c[..., None])[m]**2)**0.5)

def update_w(my_frames, w, W, b, B, P, K, minval = 1e-8, iters = 4):
    classes, rotations, pixels = W.shape
    frames = len(my_frames)
    
    # xmax[d] = (sum_tr P[d, t, r]) (sum_i K[d, i]) / sum_tr (sum_i W[t, r, i]) P[d, t, r]
    #                               (sum_i K[d, i]) / sum_tr (sum_i W[t, r, i]) P[d, t, r]
    ksums = np.sum(K[my_frames], axis=-1)
    Wsums = np.sum(W, axis=-1)
    c    = np.sum(P[my_frames] * Wsums, axis = (1,2))
    c[c<minval] = minval
    xmax = ksums / c
    xmax[xmax<minval] = minval
    
    #- gw[t, r, i] = sum_tri W[t, r, i] P[d, t, r] (K[d, i] / T[d, t, r, i] - 1)
    #    = sum_tri P[d, t, r] K[d, i] W[t, r, i] / T[d, t, r, i]  - \sum_tri W[t, r, i] P[d, t, r] 
    f  = np.empty((frames,))
    g  = np.empty((frames,))
    T  = np.empty((classes, rotations, pixels))
    PK = np.empty((classes, rotations, pixels))
    for i in range(iters) :
        f.fill(0)
        g.fill(0)
        
        for j, d in tqdm(enumerate(my_frames), total = frames, desc = f'updating fluence estimates, iteration {i}', disable = silent):
            T[:]  = w[d] + np.dot(b[d], B) / W
            PK[:] = P[d, :, :, None] * K[d]
            f[j]  = np.sum(PK / T)
            g[j] -= np.sum(PK / T**2)
        
        g[g > -minval] = -minval
        w[my_frames] += f / g * (1 - f / c)
        m = w[my_frames] > minval
        w[my_frames] = np.clip(w[my_frames], minval, xmax)    
        
        if not silent :
            print(i, np.mean((f - c)[m]**2)**0.5)


# only works for one background class
def update_b(my_frames, w, W, b, B, P, K, minval = 1e-8, iters = 4):
    # gb[d] = sum_tri P[d, t, r] K[d, i] B[l, i] / T[d, t, r, i]  - \sum_i B[l, i] 
    classes, rotations, pixels = W.shape
    frames = len(my_frames)
    
    # xmax[d, l] = sum_i K[d, i] / sum_i B[l, i]
    ksums = np.sum(K[my_frames], axis=-1)
    for t in range(b.shape[1]) :
        t_bar  = [i for i in range(b.shape[1]) if i != t]
        Bsums = np.sum(B[t], axis=-1)
        c     = Bsums
        c     = max(minval, c)
        xmax = ksums / c
        xmax[xmax<minval] = minval
        
        f  = np.empty((frames,))
        g  = np.empty((frames,))
        T  = np.empty((classes, rotations, pixels))
        PK = np.empty((classes, rotations, pixels))
        for i in range(iters) :
            f.fill(0)
            g.fill(0)
            
            for j, d in tqdm(enumerate(my_frames), total = frames, desc = f'updating background weights, iteration {i}', disable = silent):
                T[:]  = (w[d] * W  + np.dot(b[d, t_bar], B[t_bar]))/ B[t] + b[d, t] 
                PK[:] = P[d, :, :, None] * K[d]
                f[j]  = np.sum(PK / T)
                g[j] -= np.sum(PK / T**2)
            
            g[g > -minval] = -minval
            b[my_frames, t] += f / g * (1 - f / c)
            m = b[my_frames, t] > minval
            b[my_frames, t] = np.clip(b[my_frames, t], minval, xmax)    
            
            if not silent :
                print(i, np.mean((f - c)[m]**2)**0.5)

def make_W_ims(W, pixel_indices, frame_shape):
    classes, rotations, pixels = W.shape
    ims = np.zeros((classes, rotations) + frame_shape)
    for c in range(classes):
        for r in range(rotations):
            ims[c, r].ravel()[pixel_indices] = W[c, r]
    
    return ims

def make_K_ims(K, pixel_indices, frame_shape):
    ims = np.zeros((K.shape[0],) + frame_shape, dtype = np.uint8)
    for d in range(K.shape[0]):
        ims[d].ravel()[pixel_indices] = K[d]
    
    return ims
