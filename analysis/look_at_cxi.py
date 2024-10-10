import argparse
import h5py
from PyQt5 import QtGui, QtCore, QtWidgets
import pyqtgraph as pg
import numpy as np
import os
import pickle
from tqdm import tqdm
import signal

#import skimage.measure
#import skimage.segmentation

import extra_geom

from constants import PREFIX

# for (much faster) local viewing
#PREFIX='/home/andyofmelbourne/Documents/2024/p7076'
#PREFIX='/gpfs/exfel/exp/SPB/202405/p007927'
geom_fnam=f'../geom/r0120.geom'

DATA_PATH = 'entry_1/instrument_1/detector_1/data'
MASK_PATH = 'entry_1/instrument_1/detector_1/good_pixels'

def parse_cmdline_args():
    parser = argparse.ArgumentParser(description='view shots from saved hits in cxi files')
    parser.add_argument('--run', type=int, help="run number of the cxi file.")
    parser.add_argument('--cxi', type=str, help="file name of the cxi file (in saved_hits).")
    parser.add_argument('-l', '--litpixels', action='store_true', help="use litpixels to sort events.")
    parser.add_argument('-m', '--apply_mask', action='store_true', help="zero bad pixels before display")
    return parser.parse_args()

def clip_scalar(val, vmin, vmax):
    """ convenience function to avoid using np.clip for scalar values """
    return vmin if val < vmin else vmax if val > vmax else val


class Application(QtWidgets.QMainWindow):
    def __init__(self, powder, data, mask, geom, sorted_indices, litpix, long, short):
        super().__init__()
        self.Z = data.shape[0]
        self.frame_index = -999
        
        self.sorted_indices = sorted_indices
         
        self.powder = powder
        self.z_data = data
        self.z_mask = mask

        self.litpix = litpix
        self.long = long
        self.short = short
         
        self.data = np.empty(np.squeeze(data[0]).shape, dtype=np.float32)

        # get image shape
        im, centre = geom.position_modules(data[0])
        
        self.centre = centre[::-1]

        self.geom = geom
        
        self.display = np.zeros(im.shape, dtype=np.float32)
          
        self.in_replot = False
        
        self.initUI()
        

    def initUI(self):
        # Define a top-level widget to hold everything
        w = QtWidgets.QWidget()
        
        # 2D plot for the cspad and mask
        self.plot = pg.ImageView()

        # add a + at the origin
        scatter = pg.ScatterPlotItem([{'pos': self.centre, 'size': 5, 'pen': pg.mkPen('r'), 'brush': pg.mkBrush('r'), 'symbol': '+'}])
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
                
                print('sorted index', i, 'file index', j, 'litpixels', self.litpix[i], 'short diam', round(1e9 * self.short[i]),'long diam', round(1e9 * self.long[i]))
                
                if self.frame_index >= 0 :
                    self.data[:] = np.squeeze(self.z_data[j] * self.z_mask)
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
        self.geom.position_modules(self.data, out = self.display)
        if not auto :
            self.plot.setImage(self.display.T[::-1])
        else :
            self.plot.setImage(self.display.T[::-1], autoRange = False, autoLevels = False, autoHistogramRange = False)

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
        
args = parse_cmdline_args()

if args.run :
    args.cxi  = PREFIX+'/scratch/saved_hits/r%.4d_hits.cxi'%args.run
else :
    args.cxi  = f'{PREFIX}/scratch/saved_hits/{args.cxi}'


with h5py.File(args.cxi) as f:
    xyz = f['/entry_1/instrument_1/detector_1/xyz_map'][()]

    # also load powder
    powder = f['/entry_1/instrument_1/detector_1/powder'][()]


f = h5py.File(args.cxi)
data = f[DATA_PATH]

if args.apply_mask:
    print('getting mask from', MASK_PATH)
    mask = f[MASK_PATH][()]
else :
    mask = 1

if args.litpixels :
    #with h5py.File(args.cxi) as f:
    if '/entry_1/instrument_1/detector_1/hit_sigma' in f :
        litpixels = f['/entry_1/instrument_1/detector_1/hit_sigma'][()]
    elif '/entry_1/instrument_1/detector_1/hit_score' in f :
        litpixels = f['/entry_1/instrument_1/detector_1/hit_score'][()]
    else :
        litpixels = f['/entry_1/instrument_1/detector_1/lit_pixels'][()]
    sorted_indices    = np.argsort(litpixels)[::-1]
    litpix    = litpixels[sorted_indices]
    sort = True
else :
    sorted_indices  = np.arange(data.shape[0])
    litpix = np.zeros((data.shape[0],))
    sort = False


long  = f['/entry_1/sizing/long_axis_diameter'][()]
short = f['/entry_1/sizing/short_axis_diameter'][()]
long  = long[sorted_indices]
short = short[sorted_indices]


geom = extra_geom.AGIPD_1MGeometry.from_crystfel_geom(geom_fnam)

# Always start by initializing Qt (only once per application)
signal.signal(signal.SIGINT, signal.SIG_DFL) # allow Control-C
app = QtWidgets.QApplication([])
    
a = Application(powder, data, mask, geom, sorted_indices, litpix, long, short)
a.show()

## Start the Qt event loop
app.exec_()

