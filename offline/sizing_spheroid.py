import numpy as np
import scipy.constants as sc
from tqdm import tqdm

class bin_masked_array():
    """
    apply NxN binning to the last two axes of a masked array of type in[mask]
    
    Equivalent to (N=2)
        bin_mask = (mask[..., ::2, ::2]  + mask[..., 1::2, ::2]
                    mask[..., ::2, 1::2] + mask[..., 1::2, 1::2]) > 0
        
        t = in[mask]
        out = (t[..., ::2, ::2]  + t[..., 1::2, ::2]
               t[..., ::2, 1::2] + t[..., 1::2, 1::2])[bin_mask]
    """
    def __init__(self, N, mask):
        if N != 1 and N is not None :
            # bin the mask
            binned_mask  = bin(mask, N)
            binned_shape = binned_mask.shape
            binned_mask  = binned_mask > 0 
            
            binned_raveled_size = np.sum(binned_mask > 0)
            
            # get bin labels
            bin_labels_binned_raveled = np.arange(binned_raveled_size)
            
            # unravel
            bin_labels_binned = -np.ones(binned_shape, dtype = int)
            bin_labels_binned[binned_mask] = bin_labels_binned_raveled
            
            # unbin
            bin_labels = -np.ones(mask.shape, dtype = int)
            for i in range(N):
                for j in range(N):
                    bin_labels[..., i::N, j::N] = bin_labels_binned
            
            # masked bin labels
            self.bin_labels_masked = bin_labels[mask]
            self.mask = mask
        
        self.N = N
        
    def bin(self, ar):
        """
        ar = frame[mask] 
        """
        if self.N != 1 and self.N is not None :
            out = np.bincount(self.bin_labels_masked, weights = ar)
        else :
            out = ar
        return out
    
    def test(self, ar):
        # test binning: mask bin then mask
        bin1 = bin(ar * self.mask)
        temp = bin(self.mask)
        bin1 = bin1[temp>0]
        assert(np.allclose(self.bin(ar[self.mask]), bin1))
        

def bin(ar, N = 4):
    """
    bin array in 4x4 windows along last two axes
    """
    if N == 1 or N is None :
        return ar
        
    axes = (-2, -1)
    shape = list(ar.shape)
    for a in axes:
        shape[a] = shape[a]//N
    
    out = np.zeros(shape, float)
    for i in range(N):
        for j in range(N):
            out += ar[..., i::N, j::N]
    return out

def quadratic_refine(f, x, lim=3):
    N = f.shape[0]
    if N < 10 :
        return x[0], False
    
    Q = np.empty((N, 10))
    for n in range(N):
        t = x[n]
        Q[n] = [t[0]**2, t[1]**2, t[2]**2, 2*t[0]*t[1], 2*t[0]*t[2], 2*t[1]*t[2], t[0], t[1], t[2], 1]
    
    q, resid, rank, cond = np.linalg.lstsq(Q, f, rcond=None)
    
    A = np.array([[q[0], q[3], q[4]], 
                  [q[3], q[1], q[5]],
                  [q[4], q[5], q[2]]])
    b = np.array([q[6], q[7], q[8]])
    
    Ainv = np.linalg.inv(A)
    
    xmin = np.dot(Ainv, -b / 2)
     
    xlim = np.array([np.min(x, axis=0), np.max(x, axis=0)])
    #xlim = (xlim[1]-xlim[0]) * lim * xlim
    
    out = np.clip(xmin, xlim[0], xlim[1])
    #out = xmin
    
    #fit = np.dot(Q, q)

    #print(q)
    
    return out, (resid/N)**0.5/f.max()#, fit


def rotate_Rz_Rx(qs, u, v): 
    cu = np.cos(u)
    su = np.sin(u)
    cv = np.cos(v)
    sv = np.sin(v)
    out = np.zeros_like(qs)
    out[0] = qs[0] * cu - qs[1] * su * cv + qs[2] * su * sv
    out[1] = qs[0] * su + qs[1] * cu * cv - qs[2] * cu * sv
    out[2] = qs[1] * sv + qs[2] * cv 
    return out

def make_spheroid_profile(a, c, u, v, qs, scale = 1):
    """
    Thanks Filipe:
        https://www.nature.com/articles/nphoton.2014.270
    error there
    
    object = ellipsoid with semi-axes radii (a, c, a)
    F      = Fourier tranform
    
    u = out of plane rot (rotation about x-axis)
    v = in plane rot (rotation about z-axis = optical axis)
    
    F    = 4 pi a^2 c (sin H - H cos H)/H^3 
    I    = |F|^2
    H(qr) = 2 pi [ a^2 qrx^2 + c^2 qry^2 + a^2 qrz^2 ]^1/2
    
    qr = Rz(u) . Rx(v) . q
    """
    qr     = rotate_Rz_Rx(qs, v, u)
    H      = 2 * np.pi * (a**2 * qr[0]**2 + c**2 * qr[1]**2 + a**2 * qr[2]**2)**0.5
    out    = 4 * np.pi * a**2 * c * (np.sin(H) - H * np.cos(H)) / H**3
    #out[:] = out**2
    # output amplitudes instead
    return np.abs(out)

def make_q_map_solid_angle_polarisation(xyz, pixel_size, wavelength, polarisation = 'x'):
    # q = 1/wav [(x, y, z) / r - (0, 0, 1)], r = sqrt{ x^2 + y^2 + z^2 }
    r = np.sum(xyz**2, axis=0)**0.5
    q = xyz / r 
    q[2] -= 1
    q /= wavelength

    # solid angle of pixel
    # doesn't work for DSSC 
    # solid_angle = pixel_area x detector_distance / r^3
    sol = pixel_size**2 * xyz[2] / r**3
    
    # polarisation = 1 - x^2 / r^2 or 1 - y^2 / r^2
    if polarisation == 'x' :
        pol = 1 - xyz[0]**2 / r**2
    elif polarisation == 'y' :
        pol = 1 - xyz[1]**2 / r**2
    else :
        pol = 1 - (xyz[1]**2 + xyz[2]**2) / r**2 / 2
    return q, sol, pol


class Sizing():
    def __init__(self, 
            mask, xyz, pixel_size, photon_energy,      
            min_rad, max_rad, N_angles,  
            size_min, size_max, N_sizes, polarisation='x', bin_size=4):
        
        # E = hf = h c / wav
        E = photon_energy 
        wav = sc.h * sc.c / E
        print('wavelength', wav)
         
        qmap, solid, pol = make_q_map_solid_angle_polarisation(xyz, pixel_size, wav, polarisation)
        corr = solid * pol
        corr /= corr.max()
        
        r = np.sum(xyz[:2]**2, axis = 0)**0.5
        print(f'min-max radius in mask {r[mask].min():.5f} - {r[mask].max():.5f}')
        print(f'user defined min-max   {min_rad:.5f} - {max_rad:.5f}')
        
        # join pixel mask and radial range
        self.mask = mask.copy()
        self.mask[r < min_rad] = False
        self.mask[r > max_rad] = False
        
        self.pixels = np.sum(self.mask)
        print(f'number of pixels after radius masking: {self.pixels}')
        
        # apply mask to qmap
        self.qmap = qmap[:, self.mask]
        self.corr = corr[self.mask]
        
        # generate lookup 4 parameters + pixels 
        self.radii, self.dr = np.linspace(size_min/2., size_max/2, N_sizes, retstep = True)
        
        # just store the indexes
        ai = np.arange(N_sizes)
        ci = np.arange(N_sizes)
        
        # but we keep c > a to save time
        ai, ci = np.meshgrid(ai, ci, indexing='ij')
        m = ci >= ai
        ai = ai[m]
        ci = ci[m]
        ac = np.array([ai.ravel(), ci.ravel()]).T
        
        # I think we only need 0->pi/4 if we don't care about Ewald curvature
        #u = np.linspace(0, np.pi/4, N_angles)
        # what if we discard out-of-plane rotations?
        u = np.array([0.])
        #v = np.linspace(0, np.pi,   N_angles)
        v = np.pi / N_angles * np.arange(N_angles) 
        
        size = ai.size * u.size * v.size
        print(f'search space size: {size:.2e}')
        
        # make binner
        self.binner = bin_masked_array(bin_size, self.mask)

        self.pixels = self.binner.bin(self.mask[self.mask]).size
        print(f'number of pixels after radius masking and binning: {self.pixels}')
        
        # generate huge lookup table
        self.lookup = np.zeros((size, self.pixels), dtype = float)
        
        # for scaling
        self.lookup2 = np.zeros((size, self.pixels), dtype = float)

        self.errs = np.zeros((size,), dtype = float)
        
        j = 0 
        for aj, cj in tqdm(ac) :
            for uj in u :
                for vj in v :
                    p = self.make_profile(self.radii[aj], self.radii[cj], uj, vj)
                    
                    # store for later
                    self.lookup[j] = p
        
                    # scale lookup for LL
                    self.lookup2[j]  = p**2
                    self.lookup2[j] /= np.sum(self.lookup2[j])
                    
                    j += 1
        
        self.ac = ac
        self.u = u
        self.v = v
        self.wav = wav
    
    def make_profile(self, a, c, u, v, bin = True):
        # gives |F|
        p = make_spheroid_profile(a, c, u, v, self.qmap)
        
        # solid angle and polarisation correction
        p *= np.sqrt(self.corr)
        
        # bin 
        if bin :
            p = self.binner.bin(p)
        
        # normalise
        p /= np.sum(p**2)**0.5
        return p

    def quadratic_refinement(self, ai, ci, i, j, k, imin):
        # tricky to index this stuff
        
        # I want the error for a and c within +- one step
        errs = []
        acv  = []
        for x in [-1, 0, 1]:
            for y in [-1, 0, 1]:
                # find the i index where a and c are offset by x, y
                t = np.where( (self.ac == np.array([ai+x, ci+y])).all(-1) )[0]
                 
                if len(t) > 0 :
                    i = t[0]
                    for z in [-1, 0, 1]:
                        kk = (k+z) % self.v.shape[0]
                        ijk = np.ravel_multi_index((i, j, kk), (self.ac.shape[0], self.u.shape[0], self.v.shape[0]))
                        errs.append(self.errs[ijk])
                        acv.append([1e9 * self.radii[ai+x], 1e9 * self.radii[ci+y], self.v[kk]])

        # now use quadratic refinement to locate the minimum
        (a, c, v), res = quadratic_refine(np.array(errs), np.array(acv))
        a = a * 1e-9
        c = c * 1e-9
        
        if res == False :
            a = self.radii[ai]
            c = self.radii[ci]
            v = self.v[k]
        
        return a, c, v
        
    def size(self, ar, quadratic_refinement = True):
        # normalise
        t    = ar[self.mask]**0.5 / np.sum(ar[self.mask])**0.5
        
        # bin and mask
        t    = self.binner.bin(t)
        
        np.sum((self.lookup - t)**2, axis=-1, out = self.errs)
        imin = np.argmin(self.errs)
        
        i, j, k = np.unravel_index(imin, (self.ac.shape[0], self.u.shape[0], self.v.shape[0]))
        ai, ci  = self.ac[i]
        
        if quadratic_refinement :
            a, c, v = self.quadratic_refinement(ai, ci, i, j, k, imin)
        else :
            a = self.radii[ai]
            c = self.radii[ci]
            v = self.v[k]
        
        u = self.u[j]
        self.imin = imin
        
        # return average diameter
        return 2 * (2*a + c) / 3, a, c, u, v

    def size_new(self, ar, quadratic_refinement = True):
        # normalise
        #t    = ar[self.mask]**0.5 / np.sum(ar[self.mask])**0.5
        
        # bin and mask
        t    = self.binner.bin(ar[self.mask])
        
        #np.sum((self.lookup - t)**2, axis=-1, out = self.errs)
        #imin = np.argmin(self.errs)
        
        #photons = np.sum(t)

        # LL
        np.sum(-t * np.log(self.lookup2), axis=-1, out = self.errs)
        imin = np.argmin(self.errs)
        
        i, j, k = np.unravel_index(imin, (self.ac.shape[0], self.u.shape[0], self.v.shape[0]))
        ai, ci  = self.ac[i]
        
        if quadratic_refinement :
            a, c, v = self.quadratic_refinement(ai, ci, i, j, k, imin)
        else :
            a = self.radii[ai]
            c = self.radii[ci]
            v = self.v[k]
        
        u = self.u[j]
        self.imin = imin
        
        # return average diameter
        return 2 * (2*a + c) / 3, a, c, u, v


if __name__ == '__main__':
    """
    import extra_geom
    import h5py
    
    mask_file = 'r0040_good_pixels.h5'
    geom_file = 'r0035_powder.geom'
    photon_energy = 6000 * sc.e
    min_rad = 0
    max_rad = 0.02
    detector_distance = 700e-3
    size_min = 5e-9
    size_max = 50e-9
    N_sizes   = 32
    N_angles  = 32
    
    with h5py.File(mask_file) as f:
        mask = f['entry_1/good_pixels'][()]
    
    # get refined geometry
    g = extra_geom.AGIPD_1MGeometry.from_crystfel_geom(geom_file)
    xyz = g.get_pixel_positions()
    xyz = np.transpose(xyz, (3, 0, 1, 2))
    pixel_size = g.pixel_size
    xyz[2] = detector_distance

    sizing = Sizing(
                mask, xyz, pixel_size, photon_energy,      
                min_rad, max_rad, N_angles,
                size_min, size_max, N_sizes, polarisation='x', bin_size=1)
    
    with h5py.File('/home/andyofmelbourne/Documents/2024/p7076/scratch/saved_hits/r0049_hits.cxi') as f:
        frame = f['entry_1/data_1/data'][396]
    
    av_diameter, long_axis_radius, short_axis_radius, theta_x, theta_z = sizing.size(frame)
    """

    """
    f = lambda x, y: (x-0.5)**2 + (y+0.1)**2
    errs = np.zeros((3,3))
    for x in [-1, 0, 1]:
        for y in [-1, 0, 1]:
            errs[x+1, y+1] = f(x, y)
    
    print(quadratic_refine(errs))
    """
    import pickle
    xyz, errs = pickle.load(open('temp.pickle', 'rb'))
    errs = errs.ravel()
    xyz = xyz.reshape(-1, 3)

    sol, res, fit = quadratic_refine(errs, xyz)
    print(sol)
