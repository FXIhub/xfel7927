import glob
import subprocess
from constants import PREFIX
import h5py

# get all powders
fnams = glob.glob(f'{PREFIX}/scratch/powder/*.h5')

dset = '/powder'

for fnam in fnams:
    command = f'python geometry_refinement_symmetry.py {fnam} {PREFIX}/scratch/det/r0551_mask.h5 ../geom/r0600.geom -z 715e-3 -d {dset} -q -o'

    with h5py.File(fnam) as f:
        if dset in f :
            run = True
        else :
            run = False

    if run :
        print(f'running {command}')
        subprocess.run(command, shell=True, check=True, text=True)

