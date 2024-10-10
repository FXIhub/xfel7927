"""
get interesting / not interesting frames
size them then show in matplotlib image
"""
import sizing_spheroid
import h5py
import numpy as np
import extra_geom
from tqdm import tqdm

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

class Comparison_array():
    
    def __init__(self, xyz, mask, rmin, rmax, geom, wav):
        frame_shape = xyz.shape[1:]

        # mask r optional
        r = (xyz[0]**2 + xyz[1]**2)**0.5
        mask[r<rmin] = False
        mask[r>rmax] = False
         
        # x y limits
        xmin = xyz[0, mask].min()
        xmax = xyz[0, mask].max()

        ymin = xyz[1, mask].min()
        ymax = xyz[1, mask].max()

        # get image shape
        # image  = zero cropped( positioned array )
        t        = np.zeros(frame_shape, dtype = np.float32)
        t[mask]  = 1
        im, sl   = crop_to_nonzero_ND(geom.position_modules(t)[0])
        im_shape = im.shape
        self.background = np.isnan(im)
        
        # get x,y values that fill the plane (i.e. extend between panel gaps)
        x = xmin + geom.pixel_size * np.arange(im_shape[0]) 
        y = ymin + geom.pixel_size * np.arange(im_shape[1]) 
        x, y = np.meshgrid(x, y, indexing='ij')
            
        xyz_fill = np.zeros((3,) + x.shape, dtype = float)
        xyz_fill[0] = x
        xyz_fill[1] = y
        xyz_fill[2] = xyz[2].ravel()[0]
        
        # get q
        qmap, solid, pol = sizing_spheroid.make_q_map_solid_angle_polarisation(xyz_fill, geom.pixel_size, wav)
    
        self.geom = geom
        self.sl   = sl
        self.qmap = qmap
        self.corr = solid * pol
        self.mask = mask
    
    def get_image(self, ar, a, c, u, v):
        # get image of data = geometry correction then zero cropping
        data_im = self.geom.position_modules(ar.astype(float))[0][self.sl]
        
        # image of fit
        fit_full = sizing_spheroid.make_spheroid_profile(a, c, u, v, self.qmap)**2
        
        # solid angle and polarisation correction
        fit_full *= self.corr
        
        # normalise amplitude
        fit_full *= np.sum(ar[self.mask]) / np.sum(fit_full[~self.background])
          
        # now merge images
        data_im[self.background] = fit_full[self.background]
        
        # mask r optional
        #r = (x**2 + y**2)**0.5
        #data_im[r<rmin] = np.nan
        #data_im[r>rmax] = np.nan
        return data_im

# get frames 9
# fit frames
# make image arrays
# matplotlib show them in a grid

if __name__ == '__main__':
    cxi_file = '/home/andyofmelbourne/Documents/2024/p7927/scratch/saved_hits/r0437_hits.cxi'

    with h5py.File(cxi_file) as f:
        xyz  = f['/entry_1/instrument_1/detector_1/xyz_map'][()]
        mask = f['/entry_1/instrument_1/detector_1/good_pixels'][()]
        wav  = np.mean(f['/entry_1/instrument_1/source_1/photon_wavelength'][()])
    
    geom_file = '../geom/r0063.geom'
    min_rad = 0
    max_rad = 0.02
    
    geom = extra_geom.AGIPD_1MGeometry.from_crystfel_geom(geom_file)
    
    comp = Comparison_array(xyz, mask, min_rad, max_rad, geom, wav)
    
    import matplotlib.pyplot as plt
    
    with h5py.File(cxi_file) as f:
        frames = f['entry_1/data_1/data']
        a = f['/entry_1/sizing/short_axis_diameter']
        c = f['/entry_1/sizing/long_axis_diameter']
        u = f['/entry_1/sizing/theta_x']
        v = f['/entry_1/sizing/theta_z']
            
        counts = f['/entry_1/instrument_1/detector_1/photon_counts'][()]
        
        for i in np.argsort(counts)[::-1]:
            fig, ax = plt.subplots(1, 1)
            fig.set_tight_layout(True)
            fig.set_size_inches(4, 4)
    
            print(a[i], c[i], u[i], v[i])
                
            # why -v !!!???
            im = comp.get_image(frames[i], a[i]/2, c[i]/2, u[i], -v[i])
            
            ax.imshow(im**0.5, vmax = 4**0.5)
            ax.set_axis_off()
            
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.7)
            
            # place a text box in upper left in axes coords
            textstr = f"average diameter {1e9 * (a[i] + c[i])/2:.1f}nm"
            ax.text(0.5, 0.95, textstr, transform=ax.transAxes, fontsize=10,
                    horizontalalignment='center', verticalalignment='center', bbox=props)
             
            plt.show()
