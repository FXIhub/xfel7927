import matplotlib.pyplot as plt
import glob
import pickle
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


# get quad positions
qs = []
runs = []

s = 'quad_positions_r'
fnams = glob.glob(f'{s}*')

for fnam in fnams:
    i = fnam.find(s) + len(s)
    runs.append(int(fnam[i:i+4]))
    qs.append(pickle.load(open(fnam, 'rb'))) 

qs = np.array(qs) / 200e-6
runs = np.array(runs)

print(qs.shape, runs.shape)

fig, axs = plt.subplots(2, 1, sharex = True)
fig.set_size_inches(10, 10)
fig.set_tight_layout(True)

for d in range(2):
    for q in range(4):
        ax = axs[d]
        ax.scatter(runs, qs[:, q, d], s = 1, label = f'quadrant {q}')
        ax.legend()
        ax.grid(visible=True, which='both', alpha = 0.3)
        ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.yaxis.set_minor_locator(MultipleLocator(10))
        ax.set_label('run number')

axs[0].set_title('shift (pixels) refined against powder sum')
axs[0].set_ylabel('x shift (pixels)')
axs[1].set_ylabel('y shift (pixels)')

plt.show()

