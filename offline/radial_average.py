
import sys
import numpy as np
import h5py
import matplotlib.pyplot as plt


PREFIX = '/gpfs/exfel/exp/SPB/202421/p007076/scratch'
Run_id =  sys.argv[1]




with h5py.File(f'{PREFIX}/saved_hits/r00{Run_id}_hits.cxi','a') as df:
#     data = df['/entry_1/instrument_1/detector_1/data'][:]
    background = df['/entry_1/instrument_1/detector_1/background'][:]
    mask = df['/entry_1/instrument_1/detector_1/good_pixels'][:]
    xyz_map = df['/entry_1/instrument_1/detector_1/xyz_map'][:]


max_x = xyz_map[0].max()
max_y = xyz_map[1].max()

radial_bin_size = 0.01
min_rad = 0
max_rad = max(max_x,max_y)*1.5


class Radial_average():
    def __init__(self, xyz, mask, radial_bin_size, min_rad, max_rad):
        self.r = (xyz[0]**2 + xyz[1]**2)**0.5
        self.radial_bin_size = radial_bin_size
        
        self.mask = mask.copy()
        self.mask[self.r < min_rad] = False 
        self.mask[self.r > max_rad] = False 
        
        # integer label for radial bin
        self.rl = np.rint(self.r[self.mask] / radial_bin_size).astype(int).ravel()
        
        # just to speed things up a little
        self.rl = np.ascontiguousarray(self.rl)
        
        # number of pixels contributing to each bin 
        self.rcounts = np.bincount(self.rl.ravel())
        
        # for reference
#         self.rs = np.arange(np.max(self.rl)+1) * self.radial_bin_size
        
    def make_rad_av(self, ar):
        rsum = np.bincount(self.rl, ar[self.mask])
        
        # normalise, might not want to do this
        rsum /= np.clip(self.rcounts, 1, None)
        return rsum


# In[190]:


r = Radial_average(xyz_map,mask,radial_bin_size,min_rad,max_rad)


# In[191]:


ave = r.make_rad_av(background)


# In[192]:


xaxis = np.linspace(min_rad,max_rad,ave.shape[0])
ave_all_rings = np.sum(ave[5:-5])

plt.figure(figsize=(10,10))
plt.loglog(xaxis[5:],ave[5:])
plt.xlabel('radius')
plt.ylabel('background intensity')
plt.title(f'sum of all radius: {ave_all_rings}')
plt.savefig(f'{PREFIX}/background_plots/_{Run_id}_background.pdf')


with h5py.File(f'{PREFIX}/background_plots/_{Run_id}_background.h5','a') as file:
    file.create_dataset("/entry_1/instrument_1/radial_average_background",data=ave)


print('done')





