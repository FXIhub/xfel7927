from tqdm import tqdm
import numpy as np
import time
import runpy
import pathlib
import os

class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)

def load_config(path):
    p = pathlib.Path(path)
    
    print(p.absolute())
    
    # returns a dict
    config = runpy.run_path(str(p.absolute()))
    
    # convert to object (like import config)
    config = Struct(**config)
    return config
    


def plot_iter(c, r, iteration = 0):
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
    
    layout = """
        PPP012
        PPP345
        www678
        abcfff
        LLLLLL
    """
    #fig = plt.figure(constrained_layout=True)
    fig = plt.figure(tight_layout=True)
    fig.set_size_inches(20, 15)
    ax_dict = fig.subplot_mosaic(layout)
    
    # P-matrix
    ax = ax_dict["P"]
    #im = ax.imshow(r.P.T**0.5, aspect='auto', origin='lower', interpolation='nearest')
    #ax.set_ylabel('class (t)')
    #ax.set_title('P-matrix')
    # assign colour (value) to each class in subset
    print(iteration, len(r.most_likely_classes))
    frames  = np.random.randint(0, r.D, size = min(256, r.D))
    array   = np.zeros((iteration, frames.shape[0]), dtype=float)
    array[-1, :] = np.linspace(0, 1, frames.shape[0])
    for i in range(iteration-2, -1, -1):
        a = r.most_likely_classes[i][frames]
        b = r.most_likely_classes[i+1][frames]
        for j in range(array.shape[1]):
            if (a[j] - b[j]) != 0 :
                array[i, j] = np.random.random()
            else :
                array[i, j] = array[i+1, j]
    im = ax.imshow(array.T, aspect = 'auto', origin='lower', interpolation='nearest')#, cmap = 'gist_ncar_r')
    ax.set_ylabel('frame')
    ax.set_xlabel('iteration')
    ax.set_title('most likely class')
    ax.set_xticks(range(iteration))
    
    # P-plots
    """
    ax = ax_dict["l"]
    t_s = [0,1,9] 
    for t in t_s :
        ax.plot(range(r.w.shape[0]), r.P[:, t], linewidth=0.8, alpha=.5, label=f'class {t}')
    ax.set_ylabel('prob. P[t]')
    ax.legend()
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_xlabel('frame number (d)')
    """

    # w
    ax = ax_dict["w"]
    ax.plot(range(r.w.shape[0]), r.w, linewidth=0.8, alpha=1.)
    ax.set_xlabel('frame number (d)')
    ax.set_ylabel('relative fluence w[d]')
    ax.spines[['top']].set_visible(False)
    
    # b
    ax2 = ax.twinx()
    for l in range(r.b.shape[1]):
        ax2.plot(range(r.b.shape[0]), r.b[:, l], linewidth=0.8, alpha=.7, c = 'k', label=f'b class {l}')
    ax2.set_ylabel('background weight')
    ax2.spines[['top']].set_visible(False)

    # W 
    # favour = sum_d PT[t, d]
    C = r.W.shape[0]
    favour = np.sum(r.P, axis=0)
    ts     = np.argsort(favour)
    most   = ts[-3:][::-1]
    middle = ts[C//2:C//2 + 3][::-1]
    least  = ts[:3][::-1]
    
    fav     = [favour[most], favour[middle], favour[least]]
    labels  = [most, middle, least]
    classes = [r.W[most], r.W[middle], r.W[least]]
    
    image = np.empty(c.frame_shape, dtype=float)
    
    # show 3 most favoured classes and 3 least favoured class
    for i in range(3):
        for j in range(3):
            k = str(3 * i + j)
            ax = ax_dict[k]
            image.fill(0)
            image.ravel()[r.pixel_indices] = classes[i][j]
            ax.imshow(c.imshow(image)**0.2, origin='lower')
            ax.axis('off')
            ax.set_title(f'class {labels[i][j]} no. of frames {round(fav[i][j])}', fontsize=10)

    # favour bar plot
    ax = ax_dict["f"]
    ax.bar(np.arange(C), favour[ts[::-1]],  width = 1, align='edge', color='lightcoral', edgecolor='k', alpha=0.8, linewidth=1)
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_xlim([0, C])
    ax.set_title("class favour")
    ax.set_ylabel("occupancy")
    ax.set_xlabel("classes (sorted)")
    
    # show up to 3 background models
    bs = min(r.B.shape[0], 3)
    
    for l, label in enumerate(['a', 'b', 'c'][:bs]):
        ax = ax_dict[label]
        image.fill(0)
        image.ravel()[r.pixel_indices] = r.B[l]
        ax.imshow(c.imshow(image)**0.2, origin='lower')
        ax.axis('off')
        ax.set_title(f'B-class {l}')
    
    # clear remaining axes
    for l, label in enumerate(['a', 'b', 'c'][bs:]):
        ax_dict[label].set_axis_off()
    
    # expectation value plot
    ax = ax_dict["L"]
    ax.bar(np.arange(len(r.expectation_values)), r.expectation_values,  width = 1, align='edge', color='lightcoral', edgecolor='k', alpha=0.8, linewidth=1)
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_xlim([0, max(100, len(r.LL))])
    #ax.set_yscale('log')
    ax.set_title("expectation values")
    ax.set_xlabel("iterations")
    
    plt.savefig(f'{c.working_dir}/recon_{str(iteration-1).zfill(3)}.pdf')

    plt.close(fig)

def plot_all_classes(c, r, iteration = 0):
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots()
    fig.set_tight_layout(True)
    fig.set_size_inches(10, 10)
    
    # W 
    # favour = sum_d PT[t, d]
    C = r.W.shape[0]
    favour = np.sum(r.P, axis=0)
    
    image = np.empty(c.frame_shape, dtype=float)
    
    # show 3 most favoured classes and 3 least favoured class
    for i in range(C):
        image.fill(0)
        image.ravel()[r.pixel_indices] = r.W[i]
        ax.imshow(c.imshow(image)**0.2, origin='lower')
        ax.axis('off')
        ax.set_title(f'class {i} no. of frames {round(favour[i])}', fontsize=10)
         
        plt.savefig(f'{c.working_dir}/class_{str(i).zfill(3)}.pdf')

    plt.close(fig)
    
    os.system(f"pdfunite {c.working_dir}/class_*.pdf {c.working_dir}/classes.pdf")
    os.system(f"rm {c.working_dir}/class_*.pdf")
    
    

def calc_f(x, a, b, t):
    t[:] = x + b
    return np.sum(a / t)

def calc_f_df(x, a, b, t):
    t[:]  = x + b
    f  = np.sum(a / t)
    df = -np.sum(a / t**2)
    return f, df


def find_root_model(x, a, b, c, t, minval = 1e-5, tol = 1e-6, maxiters = 10, return_points = False, verb = False):
    """
    find the root of:
        f(x) - c = 0
        
        f(x) = \sum_i a[i] / (x + b[i]) 
    """
    xs = []
    
    if return_points : xs = [x]
    
    # model as u / (x + v) method
    for i in range(maxiters):
        f, g = calc_f_df(x, a, b, t)
         
        # distance to intercept
        df = f - c

        # find step size in model that produces df
        u = - f**2 / g
        v = - f / g - x
        
        xp = u / c - v
        step = xp - x
        
        x += step
        
        if return_points : xs.append(x)
        
        # check if we hit the left boundary
        if x <= minval :
            x = minval
            if verb : print('hit left boundary', minval)
            continue

        # check if we hit the right boundary
        if x <= minval :
            x = minval
            if verb : print('hit right boundary', minval)
            continue
        
        error = abs(step) / x
        if verb : print(i, 'x', x, 'step', step, 'error:', error, 'tol', tol, 'max iters', maxiters)
        
        if abs(error) < tol :
            break
    
    return x, i, xs

def find_root_Newton(x, a, b, c, t, minval = 1e-5, tol = 1e-6, maxiters = 10, return_points = False, v = False):
    """
    find the root of:
        f(x) - c = 0
        
        f(x) = \sum_i a[i] / (x + b[i]) 
    """
    xs = []
    
    if return_points : xs = [x]
    
    # newton method
    for i in range(maxiters):
        f, g = calc_f_df(x, a, b, t)
         
        # distance to intercept
        df = f - c
        
        # tangent line intersects c at:
        # f(x) = g (x - x0) + f0 = c
        # x = x0 + (c - f0) / g
        step = - df / g
         
        x += step
        
        # check if we hit the left boundary
        if x <= minval :
            x = minval
            if v : print('hit left boundary', minval)
            continue
        
        if return_points : xs.append(x)
        
        error = abs(step) / x
        if v : print(i, 'x', x, 'step', step, 'error:', error)
        
        if abs(error) < tol :
            break
    
    return x, i, xs
    
    

    
def find_root_f(x, a, b, c, minval = 1e-5, tol = 1e-6, method = "model", maxiters = 10, return_points = False, v = False):
    """
    find the root of:
        f(x) - c = 0
        
        f(x) = \sum_i a[i] / (x + b[i]) 

    xmax = 1/c \sum_i a[i]

    Warning: modifies x, a, b, c
    """
    if v: print('')
    if v: print('a-vector:', a)
    if v: print('b-vector:', b)
    if v: print('c-value :', c)
    
    # check for empty list
    if a.size == 0 :
        if v: print('empty a-vector, aborting', x)
        return x, 0, []
    
    # check for negative c (shouldn't happen)
    if c <= 0. :
        print('Warning: gradient offset <= 0!', x)
        return x, 0, []
    
    # calculate maximum value for x
    xmax = 1/c * np.sum(a)
    
    # check for zero a 
    if np.allclose(xmax, 0) :
        print('Warning: all a values = 0, aborting')
        return x, 0, []
     
    # check the right boundary
    if np.allclose(b, 0.) :
        if v : print('all b-values = 0, returning xmax', xmax)
        return xmax, 0, []
    
    # if all b[i] > 0, then it is safe to evaluate f(0)
    if np.all(b > 0.) :
        if v : print('all b-values > 0, setting minval to 0')
        minval = 0.
    
    # make a temporary array for holding x + b
    t = np.empty_like(b)
    
    # check the left boundary
    fmax = calc_f(minval, a, b, t)
    if fmax <= c :
        if v : print('c > fmax, setting x to minval', minval)
        return minval, 0, []
    
    # scale f to 1 at the mid point
    s = 1 / calc_f(xmax / 2, a, b, t)
    
    # scale a
    a[:] *= s / xmax
    
    # scale b
    b[:] /= xmax
    
    # scale c
    c *= s
    
    # scale x
    x /= xmax
    
    # make sure x is in boundary
    eps = 1e-8
    x = np.clip(x, eps, 1 - eps)
    
    # if we are far from the solution then get closer
    # get initial estimate by fitting a / (x + b) to fmax and 1
    # make sure we are to the right of the estimate
    x0 = np.clip(0.5 * (s * fmax / c - 1) / (s * fmax - 1), eps, 1-eps)
    if abs(x0 - x) > 0.3 :
        x = x0 
        if v: print('setting x to estimate:', x0)
    
    # now find the root of s f(xmax x, s a / xmax, b / xmax, s c)
    if method == "newton" :
        x, iters, xs = find_root_Newton(x, a, b, c, t, minval, tol, maxiters, return_points, v)
    
    if method == "model" :
        x, iters, xs = find_root_model(x, a, b, c, t, minval, tol, maxiters, return_points, v)
    
    if v : print('returning xmax x', xmax, x, x*xmax, 'after', iters + 1, 'iterations')
    
    return xmax * x, iters + 1, xs
    
    
    

def calculate_T_pixel(w, W, b, B, T):
    """
    T[t, d] = w[d] W[t] + sum_l b[d, l] B[l]
    """
    T[:] = w[None, :] * W[:, None] + np.sum(b * B, axis=1)
    
def calculate_T_frame(w, W, b, B, T):
    """
    T[t, i] = w W[t, i] + sum_l b[l] B[l, i]
    """
    T[:, :W.shape[1]] = w * W + np.dot(b, B)

def calculate_P(K, inds, w, W, b, B, LR, P, PT, beta, min_val = 1e-10):
    """
    LR[d, t] = sum_i (K[d, i] log(T[d,t,i]) - T[d,t,i])
    T[d,t,i] = w[d] W[d, i] + sum_l b[d, l] B[l, i]

    expectation value = sum_dt P[d, t] LR[d, t]
    """
    D = w.shape[0]
    (C, I) = W.shape
    
    expectation_value = 0.
    
    T = np.empty((C, I)) 
    for d in tqdm(range(D)):
        pixels = inds[d]
        Id     = len(inds[d])
        calculate_T_frame(w[d], W, b[d], B, T)
        T[:] = np.clip(T, min_val, None)
        LR[d]  = np.sum( K[d] * np.log(T[:, pixels]), axis=-1) 
        LR[d] -= np.sum( T, axis=-1)
        LR[d] *= beta
        
        expectation_value += np.sum(P[d] * LR[d])
    
    normalise_LR(LR, P)
    PT[:] = np.ascontiguousarray(P.T)
    return expectation_value

def normalise_LR(LR, P):
    """
    P[d,t] = exp(LR[d,t]) / sum_t exp(LR[d,t])
    """
    # normalise to produce probability matrix
    m = np.max(LR, axis=1)
    
    LR[:] = np.exp(LR - m[:, None])
    P[:]  = LR / np.sum(LR, axis=-1)[:, None]

def calculate_gradients_WB(K, w, b, T, P, gw, gb):
    """
    P must be transposed
    P[t, d] = P[class t, frame d]
    """
    # g[t, d] = (K / T - 1)
    g  = K / T - 1 
    
    # gw[t] = sum_d P[t, d] w[d] g[t, d]
    gw[:] = np.sum(P * w * g, axis = 1) 

    # gb[l] =  sum_d b[d, l] sum_t P[t, d] g[t, d]
    gb[:] = np.sum(b.T * np.sum(P * g, axis=0), axis=1)

def calculate_curvature_WB(K, w, b, T, P, gw, gb):
    c = K / T**2
    cwb = - np.sum( (np.sum(b * gb, axis=1) + w[None, :] * gw[:, None])**2 * P * c)
    return cwb

def calculate_LL_pixel(P, K, w, W, b, B, T):
    """
    calculate sum of LL terms for a single pixel
    LL = sum_{dt} P[t, d] (K[d] log T[t, d] - T[t, d])
    """
    calculate_T_pixel(w, W, b, B, T)
    #print(T.min(), T.max())
    return np.sum( P * (K * np.log(T) - T) )

class Solve_w_frame():
    """
    grad = sum_ti P[t] W[t, i] ( K[i] / T[t, i] - 1)
    
    curv = - sum_ti (grad W[t, i])^2 P[t] K[i] / T[t, i]^2
    
    T[t, i] = w W[t, i] + sum_l b[l] B[l, i]
    """

    def __init__(self, K, w, W, b, B, P, min_val = 1e-10, tol_pwk = 1e-5, tol = 1e-5, v = False, plot = False):
        self.tol     = tol
        self.min_val = min_val
        self.v       = v
        self.plot    = plot
        
        PW   = P[:, None] * W
        PWK  = PW * K
        PKW2 = PWK * W
        
        # mask by P W K 
        mask = PWK > tol_pwk
        
        self.mask = mask
        
        # precalculate gradient outside mask
        self.grad0 = - np.sum( PW[~mask] )
        
        # precalculate other stuff
        # I probably shouldn't broadcast over t so much
        self.PW   = np.ascontiguousarray(PW[mask])
        self.PWK  = np.ascontiguousarray(PWK[mask])
        self.PKW2 = np.ascontiguousarray(PKW2[mask])
        self.B    = np.ascontiguousarray((0 * W + np.sum(b[:, None] * B, axis=0))[mask])
        self.W    = np.ascontiguousarray(W[mask])
        
        self.x = w
        self.T = np.empty_like(self.W)
        
        # this is for error calculation (a bit wasteful since we don't need it)
        if v :
            self.P_full = P
            self.K_full = K
            self.B_full = np.sum(b[:, None] * B, axis=0)
            self.W_full = W
    
    def grad_curv(self, x):
        self.T[:] = x * self.W + self.B
        T = self.T
        
        g = np.sum(self.PWK / T - self.PW) + self.grad0
        
        c = - np.sum( g**2 * self.PKW2 / T**2 )
        return g, c
        
    def err(self, x):
        """
        L = sum_ti P[ti] (K[i] log T[ti] - T[ti])
        """
        T = x * self.W_full + self.B_full
        
        out  = np.sum( self.P_full[:, None] * (self.K_full * np.log(T) - T) ) 
        return out
    
    def solve(self):
        t0 = time.time()
        x = optimise_min_val_scalar(self, self.min_val, tol = self.tol, max_iters = 1000, v = self.v, plot = self.plot)
        dt = time.time() - t0
        if self.v : print('time taken: {:.2e} seconds'.format(dt))
        return x
        
class Solve_B_pixel_class():
    """
    grad = sum_td P[t, d] b[d] ( K[d] / T[t, d] - 1 )
    
    curv = - sum_td (b[d] grad)^2 P[t, d] K[d] / T[t, d]^2
    
    T[t, d] = b[d] B[l] + w[d] W[t] + sum_l' b[d, l'] B[l'] (l' != l)
    
    T[t, d], PbK, PKb2, b are big arrays to store ~4GB each
    """
    def __init__(self, K, w, W, b, B, P, l, min_val = 1e-10, tol_pbk = 1e-5, tol = 1e-5, plot = False, v = False):
        self.min_val = min_val
        self.tol = tol
        self.v = v
        self.plot = plot
        
        # mask based on P b K (same for all pixels but not classes)
        PbK = P * b[:, l] * K
        
        mask = PbK > tol

        self.mask = mask
        
        # precalculate gradient offset
        self.grad0 = - np.sum(P * b[:, l])
        
        self.x = B[l]
        
        # precalculate stuff
        self.PbK  = np.ascontiguousarray(PbK[mask])
        self.PKb2 = np.ascontiguousarray((PbK * b[:, l])[mask])
        
        # T offset
        back      = np.sum( b * B[None, :], axis = 1)
        self.T0   = np.ascontiguousarray( (w[None, :] * W[:, None] + back - b[:, l] * self.x)[mask] ) 
        self.T    = np.zeros_like(self.T0)
        
        # broadcast from d to t,d (could be expensive)
        self.b    = np.ascontiguousarray( (b[None, :, l] + 0*W[:, None])[mask])

        if v or plot: 
            self.P_full = P
            self.K_full = K
            self.b_full = b[:, l]
            self.T0_full = np.ascontiguousarray( (w[None, :] * W[:, None] + back - b[:, l] * B[l]) ) 
     
    def grad_curv(self, x):
        self.T[:] = x * self.b + self.T0
        T = self.T
        
        g = np.sum( self.PbK / T ) + self.grad0
        
        c = - np.sum( g**2 * self.PKb2 / T**2 )
        return g, c
        
    def solve(self):
        t0 = time.time()
        x = optimise_min_val_scalar(self, self.min_val, tol = self.tol, max_iters = 1000, v = self.v, plot = self.plot)
        dt = time.time() - t0
        if self.v : print('time taken: {:.2e} seconds'.format(dt))
        return x

    def err(self, x):
        """
        L = sum_td P[t, d] (K[d] log T[t, d] - T[t, d])
        """
        T = self.T0_full + x * self.b_full
        
        out  = np.sum( self.P_full * (self.K_full * np.log(T) - T) ) 
        return out
        
class Solve_b_frame():
    """
    grad = sum_ti P[t] W[t, i] ( K[i] / T[t, i] - 1)
    
    curv = - sum_ti (grad W[t, i])^2 P[t] K[i] / T[t, i]^2
    
    T[t, i] = w W[t, i] + sum_l b[l] B[l, i]
    """

    def __init__(self, K, w, W, b, B, P, min_val = 1e-10, tol_pwk = 1e-5, tol = 1e-5, v = False, plot = False):
        self.tol     = tol
        self.min_val = min_val
        self.v       = v
        self.plot    = plot
        
        PW   = P[:, None] * W
        PWK  = PW * K
        PKW2 = PWK * W
        
        # mask by P W K 
        mask = PWK > tol_pwk
        
        self.mask = mask
        
        # precalculate gradient outside mask
        self.grad0 = - np.sum( PW[~mask] )
        
        # precalculate other stuff
        # I probably shouldn't broadcast over t so much
        self.PW   = np.ascontiguousarray(PW[mask])
        self.PWK  = np.ascontiguousarray(PWK[mask])
        self.PKW2 = np.ascontiguousarray(PKW2[mask])
        self.B    = np.ascontiguousarray((0 * W + np.sum(b[:, None] * B, axis=0))[mask])
        self.W    = np.ascontiguousarray(W[mask])
        
        self.x = w
        self.T = np.empty_like(self.W)
        
        # this is for error calculation (a bit wasteful since we don't need it)
        if v :
            self.P_full = P
            self.K_full = K
            self.B_full = np.sum(b[:, None] * B, axis=0)
            self.W_full = W
    
    def grad_curv(self, x):
        self.T[:] = x * self.W + self.B
        T = self.T
        
        g = np.sum(self.PWK / T - self.PW) + self.grad0
        
        c = - np.sum( g**2 * self.PKW2 / T**2 )
        return g, c
        
    def err(self, x):
        """
        L = sum_ti P[ti] (K[i] log T[ti] - T[ti])
        """
        T = x * self.W_full + self.B_full
        
        out  = np.sum( self.P_full[:, None] * (self.K_full * np.log(T) - T) ) 
        return out
    
    def solve(self):
        t0 = time.time()
        x = optimise_min_val_scalar(self, self.min_val, tol = self.tol, max_iters = 1000, v = self.v, plot = self.plot)
        dt = time.time() - t0
        if self.v : print('time taken: {:.2e} seconds'.format(dt))
        return x
     
class Solve_W_pixel_class():
    """
    grad =   sum_d P[d] w[d] ( K[d] / T[d] - 1)
    
    curv = - sum_d (w[d] grad)^2 P[d] K[d] / T[d]^2
    
    T[d] = w[d] W + sum_l b[d, l] B[l]
    """
    
    def __init__(self, K, w, W, b, B, P, min_val = 1e-10, tol_pw = 1e-5, tol = 1e-5, plot = False):
        self.min_val = min_val
        self.tol     = tol
        self.plot    = plot
        
        # mask frames with low P w values (applies for all pixels, but not all classes)
        self.mask_Pw = (P * w) > tol_pw
        
        # precalculate gradient for frames wtih zero counts (applies for particular pixel)
        self.mask_K = np.ones( P.shape, dtype = bool )
        self.mask   = np.ones( P.shape, dtype = bool )
        
        self.update_pixel(K, w, W, b, B, P)
    
    def solve(self, v=True):
        t0 = time.time()
        x = optimise_min_val_scalar(self, self.min_val, tol = self.tol, max_iters = 1000, v = v, plot = self.plot)
        dt = time.time() - t0
        if v : print('time taken: {:.2e} seconds'.format(dt))
        return x
    
    def update_pixel(self, K, w, W, b, B, P):
        # precalculate gradient for frames wtih zero counts (applies for particular pixel)
        self.mask_K[:] = K > 0 
        
        self.grad0 = - np.sum((P*w)[self.mask_Pw & ~self.mask_K])
        
        # now combine masks
        self.mask[:] = self.mask_Pw * self.mask_K
        
        # for LL calc
        if self.plot :
            self.P2 = P[~self.mask]
            self.w2 = w[~self.mask]
            self.B2 = np.sum(b * B[None, :], axis=1)[~self.mask]

        # precalculate stuff
        self.Pw   = P[self.mask] * w[self.mask]
        self.PKw2 = P[self.mask] * w[self.mask]**2 * K[self.mask]
        self.B    = np.sum(b * B[None, :], axis=1)[self.mask]
        self.K    = K[self.mask]
        self.P    = P[self.mask]
        self.w    = w[self.mask]
        
        self.x = W
        self.T = np.zeros(self.w.shape, dtype=W.dtype)
    
    def grad_curv(self, x):
        """
        grad =   sum_d P[d] w[d] ( K[d] / T[d] - 1)
         
        curv = - sum_d (w[d] grad)^2 P[d] K[d] / T[d]^2
        
        T[d] = w[d] W[d] + sum_l b[d, l] B[l]
        """
        self.T[:] = self.w * x + self.B
        T = self.T
        
        # gradient
        g = np.sum( self.Pw * (self.K / T - 1) ) + self.grad0
        
        # curvature
        c = - np.sum( g**2 * self.PKw2 / T**2 )
        
        return g, c
    
    def err(self, x):
        """
        LL = sum_
        """
        self.T[:] = self.w * x + self.B
        out  = np.sum( self.P * (self.K * np.log(self.T) - self.T) ) 
        out -= np.sum( self.P2 * (self.w2 * x + self.B2) ) 
        return out
    
class Solve_WB_pixel():
    def __init__(self, K, w, W, b, B, P, min_val = 1e-10):
        """
        K[d], w[d], W[t], b[d, l], B[l], P[t, d]
        """
        self.K = K
        self.w = w
        self.W = W.copy() # these are temp arrays
        self.b = b
        self.B = B.copy() # these are temp arrays
        self.P = P
        self.min_val = min_val
        
        self.gw = np.empty_like(W)
        self.gb = np.empty_like(B)
        
        self.T = np.empty((W.shape[0], w.shape[0]), dtype = W.dtype)
        self.x = np.empty((W.size + B.size), dtype = W.dtype)
        self.g = np.empty((W.size + B.size), dtype = W.dtype)
        
        self.WB_to_x(W, B)
    
    def update_pixel(self, K, W, B):
        self.K[:] = K
        self.W[:] = W
        self.B[:] = B
        self.WB_to_x(W, B)

    def x_to_WB(self, x):
        self.W[:] = x[:self.W.size].reshape(self.W.shape)
        self.B[:] = x[self.W.size:].reshape(self.B.shape)
        
    def WB_to_x(self, W, B):
        self.x[:W.size] = W.ravel()
        self.x[W.size:] = B.ravel()
         
    def solve(self, v=True, plot=False):
        t0 = time.time()
        self.x = optimise_min_val(self, self.min_val, tol=1e-5, max_iters = 1000, v = v, plot = plot)
        self.x_to_WB(self.x)
        dt = time.time() - t0
        if v : print('time taken: {:.2e} seconds'.format(dt))
        return self.W, self.B

    def grad(self, x):
        self.x_to_WB(x)
        
        calculate_T_pixel(self.w, self.W, self.b, self.B, self.T)
        
        calculate_gradients_WB(self.K, self.w, self.b, self.T, self.P, self.gw, self.gb)
        self.g[:self.W.size] = self.gw
        self.g[self.W.size:] = self.gb
        
        return self.g
    
    def curv(self, d):
        self.gw[:] = d[:self.W.size].reshape(self.W.shape)
        self.gb[:] = d[self.W.size:].reshape(self.B.shape)
        
        c = calculate_curvature_WB(self.K, self.w, self.b, self.T, self.P, self.gw, self.gb)
        return c
    
    def err(self, x):
        self.x_to_WB(x)
        return calculate_LL_pixel(self.P, self.K, self.w, self.W, self.b, self.B, self.T)
        
def optimise_min_val_scalar(s, min_val, tol=1e-5, max_iters=1000, v = True, plot = False):
    x    = s.x
    mask = x > min_val
    
    for iters in range(max_iters):
        # gradient, curvature
        g, c = s.grad_curv(x)

        slope = np.sum(g**2)
        
        # if the gradient is zero and x is not at the boundary 
        # then we are done
        if (slope < tol**2) or ((g < 0) and (mask == False)) :
            if v: print('gradient within tolerance, done!')
            break
        
        # or the gradient is negative and x is at the boundary
        # then we are done
        if ((g < 0) and (mask == False)) :
            if v: print('gradient negative at boundary, done')
            break
        
        if (mask == False) and (g > 0) :    
            if v: print('positive gradient at boundary, unblocking')
            mask = True
                
        # check for non-negative curvature
        if c < 0 :
            step = - slope / c
        else :
            if v : print('Warning, non-negative curvature, taking small step in grad')
            step = 1e-5 / slope**0.5
        
        # get maximum step size before hitting min_val
        # negative gradient of unmasked x values
        if g < 0 :
            step_max = (min_val - x) / g
        else :
            step_max = np.inf
        
        if plot :
            plot_errs(x, slope, g, c, step, step_max, s.err)
        
        # take step 
        x += min(step, step_max) * g
        
        # block x if it hits boundary
        if step_max <= step:
            if v: print(f'hit lower bound, blocking')
            mask = False
            x = min_val
            
        # print error
        if v : print(iters, s.err(x), slope**0.5)
    
    return x

def optimise_min_val(s, min_val, tol=1e-5, max_iters=1000, v = True, plot = False):
    x = s.x
    mask = x > min_val
    inds = np.arange(x.size)
    
    for iters in range(max_iters):
        # gradient
        g = s.grad(x)
        g[~mask] = 0
        
        slope = np.sum(g**2)

        # if the gradient is zero we may be done
        if slope < tol**2 :
            if v: print('gradient almost zero', slope**0.5)
            
            # check if we need to release constraint
            if np.any(mask == False) :
                g = s.grad(x)
                
                i = inds[~mask][np.argmax(g[~mask])]
                if g[i] > tol :
                    if v: print('found masked value with positive gradient')
                    if v: print(f'index: {i} gradient: {g[i]}')
                    mask[i] = True
                    
                    # now continue iteration from begining
                    continue
            
            if v: print('gradient within tolerance, done!')
            break
                
        # curvature 
        c = s.curv(g)

        # check for non-negative curvature
        if c < 0 :
            step = - slope / c
        else :
            if v : print('Warning, non-negative curvature, taking small step in grad')
            step = 1e-5 / slope**0.5

        # get maximum step size before hitting min_val
        # negative gradient of unmasked x values
        m = mask & (g < 0)
        if np.any(m) :
            t = (min_val - x[m]) / g[m]
            i = np.argmin(t)
            step_max = t[i]
            step_max_index = inds[m][i]
        else :
            step_max = np.inf

        if plot :
            plot_errs(x, slope, g, c, step, step_max, s.err)
        
        # take step 
        x[mask] += min(step, step_max) * g[mask]
        
        # block x if it hits boundary
        if step_max <= step:
            if v: print(f'hit lower bound, blocking index {step_max_index}')
            mask[step_max_index] = False
            x[step_max_index] = min_val
            
        # print error
        if v : print(iters, s.err(x), slope**0.5)
    
    return x

def plot_errs(x, slope, g, c, step, step_max, err):
    import matplotlib.pyplot as plt
    alphas = step * np.linspace(0, 2, 100)
    
    # plot error profile
    err0 = err( x )
    
    alphas2 = alphas[alphas < step_max]
    errs = [err( x + a * g ) for a in alphas2]
    
    fig, ax = plt.subplots()
    ax.plot(alphas2, errs, c='k', label='numerical')
    ax.plot(alphas, err0 + alphas * slope + 0.5 * c * alphas**2, label='fit')
    
    ylim = ax.get_ylim()
    if step_max is not None :
        ax.vlines(step_max, ylim[0], ylim[1], linestyles='--', label='max_step')
    ax.vlines(step, ylim[0], ylim[1], colors='k', linestyles='--', label='step')
    ax.set_xlim([alphas.min(), alphas.max()])
    ax.legend()
    plt.show()


def update_w(P, w, W, b, B, K, inds, tol_P = 1e-5, tol = 1e-5, min_val = 1e-10, max_iters=1000):
    """
    g = sum_ti P[t] W[t, i] (K[i] / T[t, i] - 1)
      = sum_ti P[t] W[t, i] K[i] / T[t, i]  - sum_t P[t] \sum_i W[t, i]
      = sum_ti P[t] W[t, i] K[i] / (w W[t, i] + sum_l b[d, l] B[l, i])  - g0
      = sum_ti P[t] K[i] / (w + sum_l b[d, l] B[l, i] / W[t, i])  - g0
    """
    t0 = time.time()
    iter_time = 0.
    
    D = w.shape[0]
    iters_av = 0
    size_av = 0
    for d in tqdm(np.arange(D, dtype=np.int32), desc='updating fluence estimates'):
        g0   = np.dot(P[d], np.sum(W, axis=1))
        t_s  = np.where(P[d] > tol_P)[0]
        i_s  = inds[d]
        j    = np.ix_(t_s, i_s)
        a    = np.ascontiguousarray((P[d, j[0]] * K[d]).ravel())
        bb   = np.ascontiguousarray(np.tile(np.dot(b[d], B[:, i_s]), len(t_s)) / np.clip(W[j].ravel(), min_val, None))
        
        t00 = time.time()
        x, iters, xs = find_root_f(w[d], a, bb, g0, 
                                   min_val, tol, "model", max_iters)
        
        w[d] = x
        iters_av += iters
        size_av  += a.shape[0]
        iter_time += time.time() - t00
    
    total_time = time.time() - t0
    
    print('total   number of iterations:', iters_av)
    print('average number of iterations:', iters_av / D)
    print('average vector size:', size_av / D)
    print('setup time                  :', round(total_time - iter_time), 's')
    print('total time                  :', round(total_time), 's')
    print('setup time / total time     :', round(100 * (total_time - iter_time) / total_time), '%')


def update_W(PT, w, W, b, B, KT, indsT, tol_P = 1e-2, tol = 1e-5, min_val = 1e-10, max_iters=1000):
    """
    g = sum_d P[d] w[d] (K[d] / T[t] - 1)
      = sum_d P[d] w[d] K[d] / T[t] - sum_d P[d] w[d]
      = sum_d P[d] w[d] K[d] / (w[d] W + sum_l b[d, l] B[l]) - g0
      = sum_d P[d] K[d] / (W + sum_l b[d, l] B[l] / w[d]) - g0
    """
    C, I = W.shape
    D = w.shape[0]
    iters_av = 0
    size_av = 0
    iters_zero = 0
    
    t0 = time.time()
    iter_time = 0.
    
    mask0 = np.zeros((D,), dtype=bool)
    mask1 = np.zeros((D,), dtype=bool)
    for c in tqdm(np.arange(C, dtype=np.int32), desc='updating classes'):
    # testing
    #for c in tqdm([0,], desc='updating classes'):
        # find d's where P > P.max() tol and w > tol
        mask0.fill(False)
        mask0[:] = PT[c] > (PT[c].max() * tol_P)

        if c == 0 :
            verb = False
        else :
            verb = False
        
        g0   = np.dot(PT[c], w)
            
        #for i in tqdm(np.arange(I, dtype=np.int32), desc='updating classes by pixel', leave = False):
        for i in np.arange(I):
            # find indicies of indsT where P > tol
            k_s = np.where(mask0[indsT[i]])[0]
            
            # find d's where indsT and P > tol
            mask1.fill(0)
            mask1[indsT[i]] = mask0[indsT[i]]
            d_s = np.where(mask1)[0]
            
            # if the list is empty that means there are no photons with significant P 
            # for this pixel and class
            # so set to min_val 
            if len(k_s) == 0 or len(d_s) == 0 :
                W[c, i] = min_val
                continue 
            
            a  = np.ascontiguousarray(PT[c, d_s] * KT[i][k_s])
            bb = np.ascontiguousarray(np.dot(b[d_s, :], B[:, i]) / np.clip(w[d_s], min_val, None) )

            if verb :
                print('\n****************\n')
                print(f'number of photons for this class, pixel {c, i}: {KT[i][k_s]}')
                print(f'probability matrix above threshold            : {PT[c, d_s]}')
                print(f'w[d]                                          : {w[d_s]}')
                print('\n****************\n')
            
            t00 = time.time()
            x, iters, xs = find_root_f(W[c, i], a, bb, g0, 
                                       min_val, tol, "model", max_iters, v = verb)

            if False :
                print('')
                print(np.where(mask0)[0])
                print(indsT[i][k_s])
                x = np.sum(KT[i][k_s])
            
            W[c, i] = x
            iters_av += iters
            size_av  += a.shape[0]
            iter_time += time.time() - t00
            if iters == 0 :
                iters_zero += 1
    
    total_time = time.time() - t0
    
    print('total   number of iterations:', iters_av)
    print('average number of iterations:', iters_av / C / I)
    print('number of zero iterations   :', iters_zero)
    print('average vector size:', size_av / C / I)
    print('setup time                  :', round(total_time - iter_time), 's')
    print('total time                  :', round(total_time), 's')
    print('setup time / total time     :', round(100 * (total_time - iter_time) / total_time), '%')

def update_B(PT, w, W, b, B, KT, indsT, tol_P = 1e-5, tol = 1e-5, min_val = 1e-10, max_iters=1000):
    """
    g = sum_dt P[t, d] b[d] (K[d] / T[d, t] - 1)
      = sum_dt P[t, d] b[d] K[d] / T[d, t] - sum_dt P[t, d] b[d]
      = sum_dt P[t, d] b[d] K[d] / (w[d] W[t] + sum_l b[d, l] B[l]) - sum_dt P[t, d] b[d]
      = sum_dt P[t, d] b[d] K[d] / (b[d] B + w[d] W[t] + sum_l'neql b[d, l'] B[l']) - g0
      = sum_dt P[t, d] K[d] / ( B + (w[d] W[t] + sum_l'neql b[d, l'] B[l']) / b[d]) - g0
    """
    C, I = W.shape
    D = w.shape[0]
    L = B.shape[0]
    iters_av = 0
    size_av = 0
    
    iter_time = 0.

    t0 = time.time()
    
    bT = np.ascontiguousarray(b.T)
            
    PK = np.zeros((C, D), dtype = PT.dtype)
    td_mask = PK > tol_P
    for i in tqdm(np.arange(I, dtype=np.int32), desc='updating background classes by pixel'):
        # calculate P[t, d] K[d] 
        PK.fill(0)
        PK[:, indsT[i]] = PT[:, indsT[i]] * KT[i]
        
        td_mask[:] = PK > tol_P
        ts, ds = np.where(td_mask)
        
        for l in np.arange(L):
            g0   = np.sum(PT * bT[l])
             
            # find d's where indsT and P > tol
            a  = np.ascontiguousarray(PK[ts, ds])
            
            ls = [j for j in range(L) if j != l]
            back = np.dot(b[:, ls], B[ls, i])
            
            bb = np.ascontiguousarray( (w[ds] * W[ts, i] + back[ds]) / np.clip(bT[l, ds], min_val, None) )
            
            t00 = time.time()
            x, iters, xs = find_root_f(B[l, i], a, bb, g0, 
                                       min_val, tol, "model", max_iters)
            
            B[l, i] = x
            iters_av += iters
            size_av += a.shape[0]
            iter_time += time.time() - t00

    total_time = time.time() - t0
    
    print('total   number of iterations:', iters_av)
    print('average number of iterations:', iters_av / L / I)
    print('average vector size:', size_av / L / I)
    print('setup time                  :', round(total_time - iter_time), 's')
    print('total time                  :', round(total_time), 's')
    print('setup time / total time     :', round(100 * (total_time - iter_time) / total_time), '%')

def update_b(P, w, W, b, B, K, inds, tol_P = 1e-5, tol = 1e-5, min_val = 1e-10, max_iters=1000):
    """
    solve for b[d, l]
    
    g = sum_ti P[t] B[i] (K[i] / T[t, i] - 1)
      = sum_ti P[t] B[i] K[i] / (w W[t, i] + sum_l b[d, l] B[l, i]) - (sum_t P[t]) (sum_i B[i])
      = sum_ti P[t] B[i] K[i] / (b B[i] + w W[t, i] + sum_l' b[d, l'] B[l', i]) - g0
      = sum_ti P[t] K[i] / (b + w W[t, i] + sum_l' b[d, l'] B[l', i] / B[i]) - g0
    """
    D = w.shape[0]
    L, I = B.shape
    C = W.shape[0]
    
    iters_av = 0
    size_av = 0
    iter_time = 0.
    
    l_s = []
    for l in range(L):
        l_s.append([m for m in range(L) if m != l])

    t0 = time.time()
    
    for d in tqdm(range(D), desc='updating background weightings'):
        t_s  = np.where(P[d] > tol_P)[0]
        i_s  = inds[d]
        j    = np.ix_(t_s, i_s)
        a    = np.ascontiguousarray((P[d, j[0]] * K[d]).ravel())
        for l in range(L):
            g0   = np.sum(P[d]) * np.sum(B[l])
            bb   = np.ascontiguousarray( ((w[d] * W[j] + np.dot(b[d, l_s[l]], B[np.ix_(l_s[l], i_s)])) / np.clip(B[l, i_s], min_val, None) ).ravel())
            
            t00 = time.time()
            x, iters, xs = find_root_f(b[d, l], a, bb, g0, 
                                       min_val, tol, "model", max_iters)
            
            b[d, l] = x
            iters_av += iters
            size_av  += a.shape[0]
            iter_time += time.time() - t00 
    
    total_time = time.time() - t0
    
    print('total   number of iterations:', iters_av)
    print('average number of iterations:', round(iters_av / D / L, 1))
    print('average vector size         :', round(size_av / D / L))
    print('setup time                  :', round(total_time - iter_time), 's')
    print('total time                  :', round(total_time), 's')
    print('setup time / total time     :', round(100 * (total_time - iter_time) / total_time), '%')

