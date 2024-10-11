import glob
import subprocess
from constants import PREFIX
import h5py
import multiprocessing


# get all powders
fnams = glob.glob(f'{PREFIX}/scratch/powder/*.h5')

dset = '/powder'
s = 'quad_positions_r'

def calculate(fnam):
    i = fnam.find(s) + len(s)
    runs.append(int(fnam[i:i+4]))
    
    command = f'python geometry_refinement_symmetry.py {fnam} {PREFIX}/scratch/det/r0551_mask.h5 ../geom/r{run:>4}.geom -z 715e-3 -d {dset} -q -o'
    
    with h5py.File(fnam) as f:
        if dset in f :
            run = True
        else :
            run = False

    if run :
        print(f'running {command}')
        subprocess.run(command, shell=True, check=True, text=True)

if __name__ == '__main__':
    pool = multiprocessing.Pool(None)
    r = pool.map_async(calculate, fnams)
    r.wait() # Wait on the results

