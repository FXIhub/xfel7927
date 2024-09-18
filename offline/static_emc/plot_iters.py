import pickle
import numpy as np
import utils
import glob
import os
import sys

from pypdf import PdfWriter
import glob

def pdfunite(pattern, output, remove=False):
    fnams = sorted(glob.glob(pattern))
    merger = PdfWriter()
    for fnam in fnams:
        merger.append(fnam)
    merger.write(output)
    merger.close()

    if remove:
        os.system(f'rm {pattern}')
    

c = utils.load_config(sys.argv[1])

class Empty():
    pass

fnams = np.sort(glob.glob(c.working_dir + '/recon_*.pickle'))

for fnam in fnams:
    a = pickle.load(open(fnam, 'rb'))
    
    utils.plot_iter(c, a, a.iterations)

#os.system(f"pdfunite {c.working_dir}/recon_*.pdf {c.working_dir}/recon.pdf")
utils.plot_all_classes(c, a, a.iterations)

pdfunite(f'{c.working_dir}/recon_*.pdf', f'{c.working_dir}/recon.pdf', remove = True)
pdfunite(f'{c.working_dir}/class_*.pdf', f'{c.working_dir}/classes.pdf', remove = True)

