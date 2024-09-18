import argparse
import h5py
from PyQt5 import QtGui, QtCore, QtWidgets
import pyqtgraph as pg
import numpy as np
import os
import pickle
from tqdm import tqdm
import signal

import skimage.measure
import skimage.segmentation

#from constants import PREFIX

# for (much faster) local viewing
PREFIX='/home/andyofmelbourne/Documents/2024/p7076'

DATA_PATH = 'entry_1/instrument_1/detector_1/data'
MASK_PATH = 'entry_1/instrument_1/detector_1/mask'

def parse_cmdline_args():
    parser = argparse.ArgumentParser(description='view shots from saved hits in cxi files')
    parser.add_argument('cxi', type=int, help="run number of the cxi file.")
    parser.add_argument('-l', '--litpixels', action='store_true', help="use litpixels to sort events.")
    parser.add_argument('-m', '--apply_mask', action='store_true', help="zero bad pixels before display")
    return parser.parse_args()

def clip_scalar(val, vmin, vmax):
    """ convenience function to avoid using np.clip for scalar values """
    return vmin if val < vmin else vmax if val > vmax else val

class Application(QtWidgets.QMainWindow):
    def __init__(self, powder, data, mask, sorted_indices, indices_image_space, background_mask, xy_map, im_shape, d, litpix):
        super().__init__()
        self.Z = data.shape[0]
        self.frame_index = -999
        
        self.sorted_indices = sorted_indices
         
        self.xmin, self.ymin, self.dx = d
        
        self.im_shape = im_shape

        self.powder = powder
        self.z_data = data
        self.z_mask = data

        self.litpix = litpix
         
        self.data = np.empty(np.squeeze(data[0]).shape, dtype=np.float32)
        
        self.pixel_map       = indices_image_space
        self.background_mask = background_mask

        self.display = np.zeros(background_mask.shape, dtype=np.float32)
          
        self.in_replot = False
        
        self.initUI()
        

    def initUI(self):
        # Define a top-level widget to hold everything
        w = QtWidgets.QWidget()
        
        # 2D plot for the cspad and mask
        self.plot = pg.ImageView()

        # add a + at the origin
        # x=0, i = -xmin / dx + 0.5
        i0 = -self.xmin / self.dx + 0.5
        j0 = -self.ymin / self.dx + 0.5
        scatter = pg.ScatterPlotItem([{'pos': (i0, j0), 'size': 5, 'pen': pg.mkPen('r'), 'brush': pg.mkBrush('r'), 'symbol': '+'}])
        self.plot.addItem(scatter)
        
        if self.Z > 1 :
            # add a z-slider for image selection
            z_sliderW = pg.PlotWidget()
            z_sliderW.plot(self.litpix, pen=(255, 150, 150))
            z_sliderW.setFixedHeight(100)
            
            # vline
            self.bounds = [-1, self.Z-1]
            self.vline = z_sliderW.addLine(x = 0, movable=True, bounds = self.bounds)
            self.vline.setValue(-1) 
            self.vline.sigPositionChanged.connect(self.replot_frame)
        
        vbox = QtWidgets.QVBoxLayout()
        vbox.addWidget(self.plot)
        vbox.addWidget(z_sliderW)
    
        ## Display the widget as a new window
        w.setLayout(vbox)
        self.setCentralWidget(w)
        self.resize(800, 480)
        
        self.replot_frame(True)
        

    def replot_frame(self, auto=False):
        if self.in_replot:
            return
        try:
            self.in_replot = True
            i = int(self.vline.value())
            if self.frame_index != i :
                self.frame_index = i
                j = self.sorted_indices[self.frame_index]
                
                print(i, j, self.litpix[i], self.data.dtype)
                
                if self.frame_index >= 0 :
                    self.data[:] = np.squeeze(self.z_data[j] * self.z_mask[j])
                elif self.frame_index == -1 :
                    self.data[:] = np.squeeze(self.powder)
                
                self.updateDisplayRGB(auto)
        finally:
            self.in_replot = False

    def updateDisplayRGB(self, auto = False):
        """
        Make an RGB image (N, M, 3) (pyqt will interprate this as RGB automatically)
        with masked pixels shown in blue at the maximum value of the cspad. 
        This ensures that the masked pixels are shown at full brightness.
        """
        self.display[~self.background_mask] = self.data.ravel()[self.pixel_map]
        if not auto :
            self.plot.setImage(self.display.reshape(self.im_shape))
        else :
            self.plot.setImage(self.display.reshape(self.im_shape), autoRange = False, autoLevels = False, autoHistogramRange = False)

    def keyPressEvent(self, event):
        super(Application, self).keyPressEvent(event)
        key = event.key()
        
        if key == QtCore.Qt.Key_Left :
            ind = clip_scalar(self.frame_index - 1, self.bounds[0], self.bounds[1]-1)
            self.vline.setValue(ind)
            self.replot_frame()
        
        elif key == QtCore.Qt.Key_Right :
            ind = clip_scalar(self.frame_index + 1, self.bounds[0], self.bounds[1]-1)
            self.vline.setValue(ind)
            self.replot_frame()
        
    

def generate_pixel_lookup(xyz, oversampling = 1):
    # choose xy bounds
    xmin = xyz[0].min()
    xmax = xyz[0].max()
    
    ymin = xyz[1].min()
    ymax = xyz[1].max()
    
    # choose sampling
    dx = 177e-6 / oversampling
    
    shape = (int( (xmax-xmin)/dx ) + 2, int( (ymax-ymin)/dx ) + 2)

    # pixel coordinates in im
    ss = np.round((xyz[0] - xmin) / dx).astype(int)
    fs = np.round((xyz[1] - ymin) / dx).astype(int)
    
    i, j = np.indices(shape)
    xy_map = np.empty((2,) + shape, dtype=float)
    xy_map[0] = dx * i
    xy_map[1] = dx * j
    
    # now use pixel indices as labels
    i = np.arange(xyz[0].size)
     
    # create an image of the data raveled indices
    im = -np.ones(shape, dtype=int)
    im[ss.ravel(), fs.ravel()] = i
    
    # label image
    # problem, the labells dont equal i
    #l = skimage.measure.label(im, background=-1)
    
    # expand by oversampling rate (to fill gaps)
    l = skimage.segmentation.expand_labels(im+1, distance = oversampling)
        
    # set background mask
    background_mask = (l.ravel()==0).copy()
    
    # now subtract 1 from labels to turn them into pixel indices
    l -= 1
    
    indices_image_space = l.ravel()[~background_mask].copy()
    
    # now to map data to 2D image we have:
    # im[~background] = data.ravel()[indices_image_space]
    return indices_image_space, background_mask, shape, (xmin, ymin, dx), im, i, l

args = parse_cmdline_args()
args.run  = args.cxi
args.cxi  = PREFIX+'/scratch/saved_hits/r%.4d_hits.cxi'%args.cxi

if args.litpixels :
    with h5py.File(args.cxi) as f:
        litpixels = f['/entry_1/instrument_1/detector_1/lit_pixels'][()]
        sorted_indices    = np.argsort(litpixels)[::-1]
        litpix    = litpixels[sorted_indices]
        sort = True
else :
    sort = False

with h5py.File(args.cxi) as f:
    xyz = f['/entry_1/instrument_1/detector_1/xyz_map'][()]

    # also load powder
    powder = f['/entry_1/instrument_1/detector_1/background'][()]


f = h5py.File(args.cxi)
data = f[DATA_PATH]

if args.apply_mask:
    mask = f[MASK_PATH]
else :
    mask = np.ones((data.shape[0],), dtype = int)

# generate pixel lookup
indices_image_space, background_mask, im_shape, (xmin, ymin, dx), im, i, l = generate_pixel_lookup(xyz)

# Always start by initializing Qt (only once per application)
signal.signal(signal.SIGINT, signal.SIG_DFL) # allow Control-C
app = QtWidgets.QApplication([])
    
a = Application(powder, data, mask, sorted_indices, indices_image_space, background_mask, xyz, im_shape, (xmin, ymin, dx), litpix)
a.show()

## Start the Qt event loop
app.exec_()

