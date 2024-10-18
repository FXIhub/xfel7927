# read in photons from each h5 file
# sort by photon count
# display

# global index is seqential list
# index events by global_index -> fnam, file_index pairs
import h5py
import numpy as np
import pyqtgraph as pg
import extra_geom
from tqdm import tqdm
from PyQt5 import QtGui, QtCore, QtWidgets
from pyqtgraph.graphicsItems.InfiniteLine import InfiniteLine
import signal
import os
import sys
import h5py
import extra_geom
import math
from collections import defaultdict

def update_h5(f, key, value, compression = False):
    if key in f and f[key].shape != value.shape:
        del f[key]
        
    if key in f and f[key].shape == value.shape:
        f[key][...] = value
    
    if key not in f:
        if compression :
            f.create_dataset(key, data = value, compression='gzip', compression_opts=1, shuffle=True, chunks=True)
        else :
            f.create_dataset(key, data = value)

# get frames from cxi file
# get models from recon file
# show models as tiled 2D image
# click on models to show frames

PREFIX = os.environ['EXP_PREFIX']

cxi_fnam   = f'{PREFIX}/scratch/saved_hits/Ery_size_filtered.cxi'
recon_fnam = f'{PREFIX}/scratch/2D-EMC/Ery/recon.h5'
geom_fnam  = '../../geom/r0600.geom'

#cxi_fnam   = '/home/andyofmelbourne/Documents/2024/p7927/scratch/saved_hits/Ery_size_filtered.cxi'
#recon_fnam = '/home/andyofmelbourne/Documents/2024/p7927/scratch/2D-EMC/Ery/recon.h5'
#geom_fnam  = '../../geom/r0600.geom'

geom = extra_geom.AGIPD_1MGeometry.from_crystfel_geom(geom_fnam)
with h5py.File(cxi_fnam) as f:
    frame = f['entry_1/data_1/data'][0]
    im, centre = geom.position_modules(frame)

frame_shape = frame.shape
frame_dtype = frame.dtype
image_shape = im.shape
image_dtype = np.uint8


# indices for each class
class_indices = {}
with h5py.File(recon_fnam) as f:
    most_likely_class = f['/iterations/most_likely_class'][-1].astype(int)
    models      = f['models'][()]

class_shape = models.shape[1:]
classes     = models.shape[0]

occupancy = np.bincount(most_likely_class)
for c in range(classes):
    class_indices[c] = np.where(most_likely_class == c)[0]

# make class image
n = math.ceil(classes**0.5)
classes_im = np.zeros((n * models.shape[1], n * models.shape[0]), dtype = np.float32)

N, M = class_shape
for i in range(n):
    for j in range(n):
        c = n*i + j
        if c < classes :
            classes_im[i*N: (i+1)*N, j*M: (j+1)*M] = models[c]

def load_class_images(c, ims = None):
    if ims is None :
        ims = np.zeros((len(class_indices[c]),) + image_shape, dtype = image_dtype)
    with h5py.File(cxi_fnam) as f:
        data = f['entry_1/data_1/data']
        for i, d in tqdm(enumerate(class_indices[c]), total = len(class_indices[c]), desc='loading frames and applying geometry'):
            frame[:] = data[d]
            ims[i]   = geom.position_modules(frame)[0]

            if i == (ims.shape[0]-1) :
                break
    return ims, len(class_indices[c])
    

# Get key mappings from Qt namespace
qt_keys = (
    (getattr(QtCore.Qt, attr), attr[4:])
    for attr in dir(QtCore.Qt)
    if attr.startswith("Key_")
)
keys_mapping = defaultdict(lambda: "unknown", qt_keys)


class ImageView(pg.ImageView):
    def __init__(self, *args, **kwargs):
        super(ImageView, self).__init__(*args, **kwargs)
        
        print('press "f" to display class images of last selection')
        print('press "s" to save selection to cxi file and good_classes.pickle')
        
        self.ims = np.zeros((2000,) + image_shape, dtype = image_dtype)
        
        self.last_selected = None
        self.selection = np.zeros((classes,), dtype = bool)
        self.pos       = np.zeros((classes,2), dtype = float)
        c = np.arange(classes)
        self.pos[:, 0] = M * (c %  n)   + M/2 + 0.5
        self.pos[:, 1] = N * (c // n) + N/2 + 0.5
        
        # set hover: fill with grey
        self.scatter_hover = pg.ScatterPlotItem(size=64, pen=None, brush=None, symbol='s', pxMode=False, hoverBrush = pg.mkBrush(255, 255, 255, 100), hoverable=True)
        #spots = [{'pos': [M * (c % n) + M/2 + 0.5, N * (c//n) + N/2 + 0.5], 'data': c} for c in range(classes)]
        #self.scatter_hover.addPoints(spots)
        self.addItem(self.scatter_hover)
            
        self.update_selection()
        
        self.scatter_hover.sigClicked.connect(self.clicked)
    
    def update_selection(self):
        spots = []
        for c in range(classes):
            pen  = pg.mkPen('g') if self.selection[c] else None
            spot = {'pos': self.pos[c], 'pen': pen, 'data': c}
            spots.append(spot)
        
        self.scatter_hover.setData(spots)
                
    def clicked(self, points, ev):
        for p in ev:
            c = p.data()
            self.selection[c]  = ~self.selection[c]
            if self.selection[c] : self.last_selected = c
            self.update_selection()
            
        
    def keyPressEvent(self, event):
        super(ImageView, self).keyPressEvent(event)
        key = keys_mapping[event.key()]
        print("key press", key)
        
        if key == 'F' and self.last_selected is not None:
            self.ims, D = load_class_images(self.last_selected, self.ims)
            pg.show(self.ims[:D])

        if key == 'S' and np.any(self.selection):
            print(f'saving selection to {cxi_fnam}')
            good_classes = np.where(self.selection)[0]
            good_frames  = np.concatenate([class_indices[c] for c in good_classes])
            with h5py.File(cxi_fnam) as f:
                Nframes = f['entry_1/data_1/data'].shape[0]
            is_sample_hit = np.zeros((Nframes,), dtype = bool)
            is_sample_hit[good_frames] = True
            
            out = {'good_classes': good_classes, 'good_frames': good_frames, 'is_sample_hit': is_sample_hit, 'class': most_likely_class}
            
            with h5py.File(cxi_fnam, 'r+') as f:
                for k, v in out.items(): 
                    key = f'/manual_class_selection/{k}'
                    update_h5(f, key, v, compression = True)
                
            # write to pickle as well
            import pickle 
            fnam = 'class_selection.pickle'
            print(f'saving selection to {fnam}')
            pickle.dump(out, open(fnam, 'wb'))

            print('total number of selected frames:', good_frames.shape[0])
                
        


signal.signal(signal.SIGINT, signal.SIG_DFL) # allow Control-C

# Always start by initializing Qt (only once per application)
app = QtWidgets.QApplication([])

imageItem  = pg.ImageItem(classes_im**0.2, axisOrder='row-major')

plot = ImageView(imageItem = imageItem)
plot.show()


## Start the Qt event loop
app.exec_()
