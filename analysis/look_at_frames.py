import argparse
import h5py
from PyQt5 import QtGui, QtCore, QtWidgets
import pyqtgraph as pg
import numpy as np
import os
import pickle
from tqdm import tqdm
import signal

import extra_geom

from constants import PREFIX, VDS_DATASET 

# for (much faster) local viewing
geom_fnam=f'../geom/r0600.geom'

def parse_cmdline_args():
    parser = argparse.ArgumentParser(description='view frames from saved VDS files')
    parser.add_argument('run', type=int, help="run number of the VDS file.")
    parser.add_argument('-l', '--litpixels',  action='store_true', help="use litpixels to sort events.")
    parser.add_argument('-p', '--photon_counts',  action='store_true', help="use total photon counts to sort events.")
    parser.add_argument('-m', '--count_in_mask',  action='store_true', help="use litpixels / photon_counts inside hit mask to sort events.")
    args = parser.parse_args()
    args.vds_file    = PREFIX+'scratch/vds/r%.4d.cxi' % args.run
    args.events_file = PREFIX+'scratch/events/' + os.path.splitext(os.path.basename(args.vds_file))[0] + '_events.h5'
    return args

def clip_scalar(val, vmin, vmax):
    """ convenience function to avoid using np.clip for scalar values """
    return vmin if val < vmin else vmax if val > vmax else val


class Application(QtWidgets.QMainWindow):
    def __init__(self, data, geom, sorted_indices, sort_data):
        super().__init__()
        self.Z = data.shape[0]
        
        self.data = data
        
        self.sorted_indices = sorted_indices
        self.sort_data      = sort_data
         
        # get image shape
        im, centre = geom.position_modules(data[0])
        
        self.centre = centre[::-1]
        
        self.geom = geom
        
        self.display = np.zeros(im.shape, dtype=np.float32)
        self.display[:] = np.nan
          
        self.in_replot = False
        self.frame_index = -1
        self.file_index  = -1
        
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
            z_sliderW.plot(self.sort_data, pen=(255, 150, 150))
            z_sliderW.setFixedHeight(100)
            
            # vline
            self.bounds = [0, self.Z-1]
            self.vline = z_sliderW.addLine(x = 0, movable=True, bounds = self.bounds)
            self.vline.setValue(0) 
            self.vline.sigPositionChanged.connect(self.replot_frame)
        
        vbox = QtWidgets.QVBoxLayout()
        vbox.addWidget(self.plot)
        vbox.addWidget(z_sliderW)

        self.replot_frame(True)
    
        ## Display the widget as a new window
        w.setLayout(vbox)
        self.setCentralWidget(w)
        self.resize(800, 480)
        

    def replot_frame(self, auto=False):
        if self.in_replot:
            return
        try:
            self.in_replot = True
            i = int(self.vline.value())
            if self.frame_index != i :
                self.frame_index = i
                self.file_index  = self.sorted_indices[self.frame_index]
                
                print()
                print('sorted index: ', self.frame_index)
                print('file index  : ', self.file_index)
                print('sort data   : ', self.sort_data[self.frame_index])
                
                self.updateDisplayRGB(auto)
        finally:
            self.in_replot = False

    def updateDisplayRGB(self, auto = False):
        """
        Make an RGB image (N, M, 3) (pyqt will interprate this as RGB automatically)
        with masked pixels shown in blue at the maximum value of the cspad. 
        This ensures that the masked pixels are shown at full brightness.
        """
        self.geom.position_modules(self.data[self.file_index], out = self.display)
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

if args.litpixels :
    if args.count_in_mask :
        sort_dataset = 'litpixels_mask'
    else :
        sort_dataset = 'litpixels'
elif args.photon_counts :
    if args.count_in_mask :
        sort_dataset = 'total_intens_mask'
    else :
        sort_dataset = 'total_intens'
else :
    sort_dataset = None
    
print(f'sorting frames by dataset {sort_dataset} in {args.events_file}')

with h5py.File(args.events_file) as f:
    if sort_dataset :
        sort_data      = f[sort_dataset][()]
        sorted_indices = np.argsort(sort_data)[::-1]
        sort_data      = sort_data[sorted_indices]
    else :
        sort_data      = np.arange(f[sort_dataset].shape[0])
        sorted_indices = np.arange(f[sort_dataset].shape[0])

geom = extra_geom.AGIPD_1MGeometry.from_crystfel_geom(geom_fnam)

data = h5py.File(args.vds_file, 'r')[VDS_DATASET]

# Always start by initializing Qt (only once per application)
signal.signal(signal.SIGINT, signal.SIG_DFL) # allow Control-C
app = QtWidgets.QApplication([])

pg.setConfigOption('background', pg.mkColor(0.1))
pg.setConfigOption('foreground', 'w')
pg.setConfigOptions(antialias=True)
    
a = Application(data, geom, sorted_indices, sort_data)
a.show()

## Start the Qt event loop
app.exec_()

