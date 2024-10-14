import h5py
import numpy as np
import scipy.constants as sc
from tqdm import tqdm
import os
from emc_2d import utils
import sys
import shutil

np.random.seed(1)

"""
(base) :~/Documents/git_repos/2D-EMC$ h5ls -r ~/Documents/2023/P3004-take-2/gold/hits_r0087.cxi 
/                        Group
/entry_1                 Group
/entry_1/cellId          Dataset {59087}
/entry_1/data_1          Group
/entry_1/data_1/data     Soft Link {/entry_1/instrument_1/detector_1/data}
/entry_1/experiment_identifier Dataset {59087}
/entry_1/instrument_1    Group
/entry_1/instrument_1/data_1 Group
/entry_1/instrument_1/detector_1 Group
/entry_1/instrument_1/detector_1/data Dataset {59087, 16, 128, 512}
/entry_1/instrument_1/detector_1/good_pixels Dataset {16, 128, 512}
/entry_1/instrument_1/detector_1/xyz_map Dataset {3, 16, 128, 512}
/entry_1/instrument_1/name Dataset {SCALAR}
/entry_1/pulseId         Dataset {59087}
/entry_1/sample_1        Group
/entry_1/sample_1/name   Dataset {SCALAR}
/entry_1/trainId         Dataset {59087}
/misc                    Group
/static_emc              Group
/static_emc/class        Dataset {59087}
/static_emc/good_classes Dataset {219}
/static_emc/good_hit     Dataset {59087}
"""

# load configuration file
config = utils.load_config(sys.argv[1])

# make output directory
if len(sys.argv) == 3 :
    out = sys.argv[2]
    write_output = True
    if not os.path.exists(out) :
        print('creating reconstruction directory (can be renamed later)', out)
        os.mkdir(out)
else :
    write_output = True
    for n in range(1, 1000):
        out = 'recon_%.4d'%n
        if not os.path.exists(out) :
            print('creating reconstruction directory (can be renamed later)', out)
            os.mkdir(out)
            break



# make data file
# --------------
output = out + '/data.cxi'

# returns everything except frames
#data = utils.write_data_file(output, config['data'], config['model_length'], config['sampling'], config['max_frames'], config['frame_slice'], config['filter_by'], config['filter_value'])


# add some datasets to the cxi file for this script to work
for fnam in config['data'] :
    with h5py.File(fnam, 'r+') as f:
        key1 = '/entry_1/instrument_1/detector_1/good_pixels'
        key2 = '/entry_1/instrument_1/detector_1/mask'
        if key2 not in f :
            f[key2] = f[key1][()]

        key1 = '/entry_1/instrument_1/detector_1/xyz_map'
        key2 = '/entry_1/instrument_1/detector_1/distance'
        if key2 not in f :
            f[key2] = np.mean(f[key1][2])

        key1 = '/entry_1/instrument_1/source_1/pulse_energy'
        key2 = '/entry_1/instrument_1/source_1/energy'
        if key2 not in f :
            f[key2] = np.mean(f[key1][()])



# determine rmax and thus which pixels to select from config['data']
with h5py.File(config['data'][0], 'r') as f:
    mask = f['/entry_1/instrument_1/detector_1/mask'][()]
    xyz  = f['/entry_1/instrument_1/detector_1/xyz_map'][()]
    dx   = f['/entry_1/instrument_1/detector_1/x_pixel_size'][()]

# hack to test effecct of geometry on classes
hack = True
if hack :
    # offset from average of powder refinements
    quads = np.array([[ -6.8, -10. ],
           [ -7.1,   2.7],
           [ -8.3,   3.9],
           [  2.7,  -3.9]])
    quads *= 200e-6 

    xyz0 = xyz.copy()
    for q in range(quads.shape[0]):
        xyz[0, 4*q: 4*(q+1)] += quads[q, 0]
        xyz[1, 4*q: 4*(q+1)] += quads[q, 1]
    xyz[2] = 715e-3

    # check
    for d in range(0, 1, 2):
        for m in range(16):
            print(f'dimension {d} module {m} mean change (pixels) {np.mean(xyz[d, m] - xyz0[d, m])/200e-6}')
    
# determine model voxel locations
dx_model = dx * config['sampling']

# location of zero pixel
i0 = np.float32(config['model_length'] // 2)

x  = dx_model * (np.arange(config['model_length']) - i0)

points_I = (x.copy(), x.copy())

if 'rmax' in config and config['rmax']:
    rmax = config['rmax']
else :
    rmax = np.abs(x).max() - 2 * dx

r = (xyz[0]**2 + xyz[1]**2)**0.5
mask[r >= rmax] = False


# write meta data
with h5py.File(config['data'][0], 'r') as f:
    data_dtype  = f['entry_1/data_1/data'].dtype
    frame_shape = f['entry_1/data_1/data'].shape[1:]
    
    datasets = ['/entry_1/instrument_1/detector_1/mask',
                '/entry_1/instrument_1/detector_1/xyz_map',
                '/entry_1/instrument_1/detector_1/x_pixel_size',
                '/entry_1/instrument_1/detector_1/y_pixel_size',
                '/entry_1/instrument_1/detector_1/pixel_area',
                '/entry_1/instrument_1/detector_1/distance',
                '/entry_1/instrument_1/detector_1/background',
                '/entry_1/instrument_1/source_1/energy']
    
    file_metadata = {}
    for d in datasets:
        v = f[d][()]

        # override xyz_map
        if hack and 'xyz_map' in d:
            v = xyz.copy()
        
        # assume any dataset with frame shape as the last dimensions
        # need to be masked (including mask)
        if v.shape[-len(frame_shape):] == frame_shape :
            v = v[..., mask]
        
        file_metadata[d] = v.copy()
    
    file_metadata['/entry_1/instrument_1/detector_1/pixel_indices'] = np.where(mask.ravel())[0]
    file_metadata['/entry_1/instrument_1/detector_1/frame_shape'] = frame_shape
    
    if write_output :
        with h5py.File(output, 'w') as g:
            for d, v in file_metadata.items():
                print('writing dataset', d,'from first input file')
                if d in g :
                    g[d][...] = v
                else :
                    g[d] = v

# get location of good frames in datasets
frames = {}
number_of_frames = 0 
for fnam in config['data'] :
    with h5py.File(fnam, 'r') as f:
        if config['filter_by'] : 
            frames[fnam] = np.where(f[config['filter_by']][()] == config['filter_value'])[0]
        else :
            frames[fnam] = np.arange(f['/entry_1/experiment_identifier'].shape[0])
        
        number_of_frames += len(frames[fnam])
        
        if number_of_frames >= config['max_frames'] :
            frames[fnam] = frames[fnam][:config['max_frames'] - number_of_frames]
            number_of_frames = config['max_frames']
            break

print('found', number_of_frames, 'to write to output file')
file_metadata['frames'] = number_of_frames

pixels = out_shape = (np.sum(mask),)

# write sparse photons (data) to file (dimension = 1)
# write photon locations (inds) to file (dimension = 1)
# write number of photons (photons) to file (dimension = number of frames)
if write_output :
    frame_index = 0
    photon_index = 0
    with h5py.File(output, 'a') as g:
        # event identifiers
        cellids  = []
        pulseids = []
        vdsids   = []
        fnams    = []
        file_index = []
        
        # assume all data fits in memory
        # writing and resizing on individual frames is too expensive (memory + time)
        inds_fnam = []
        data_fnam = []
        photons_fnam = []
        litpix_fnam = []
        for i, fnam in enumerate(frames) :
            fnams.append(fnam)
            
            with h5py.File(fnam, 'r') as f:
                data_file     = f['entry_1/data_1/data']
                cellids_file  = f['/entry_1/cellId'][()]
                pulseids_file = f['/entry_1/pulseId'][()]
                vdsids_file   = f['/entry_1/experiment_identifier'][()]
            
                for d in tqdm(frames[fnam], desc=f'loading data: {fnam}'):
                    frame = data_file[d][mask]
                    inds_fnam.append(np.where(frame>0)[0])
                    data_fnam.append(frame[inds_fnam[-1]])
                    photons_fnam.append(np.sum(data_fnam[-1]))
                    litpix_fnam.append(len(inds_fnam[-1]))

                    cellids.append(cellids_file[d])
                    pulseids.append(pulseids_file[d])
                    vdsids.append(vdsids_file[d])
                    file_index.append(i)
         
        litpix_fnam  = np.array(litpix_fnam)
        inds_fnam    = np.concatenate(inds_fnam)
        data_fnam    = np.concatenate(data_fnam)
        photons_fnam = np.array(photons_fnam)

        cellids    = np.array(cellids)
        pulseids   = np.array(pulseids)
        vdsids     = np.array(vdsids)
        file_index = np.array(file_index)
        
        dt = h5py.special_dtype(vlen=str)
        fnams = np.array(fnams, dtype=dt)

        print()
        print(f'writing sparse photons to {output}')
        g['entry_1/data_1/data'] = data_fnam
        g['entry_1/data_1/inds'] = inds_fnam.astype(np.int32)
        g['entry_1/data_1/photons'] = photons_fnam.astype(np.int32)
        g['entry_1/data_1/litpix'] = litpix_fnam.astype(np.int32)

        g['entry_1/cellId'] = cellids
        g['entry_1/pulseId'] = pulseids
        g['entry_1/vds_frame_index'] = vdsids
        g['entry_1/vds_file_index'] = file_index
        g['entry_1/vds_file_names'] = fnams

J         = config['model_length']
classes   = config['classes']
rotations = config['rotations']
dtype     = config['dtype']

# make reconstruction flie
# ------------------------
recon_file = out + '/recon.h5'

# solid angle and polarisation correction
C = utils.solid_angle_polarisation_factor(file_metadata['/entry_1/instrument_1/detector_1/xyz_map'], 
                                          file_metadata['/entry_1/instrument_1/detector_1/pixel_area'], 
                                          polarisation_axis = config['polarisation_axis'])


# rotation matrices
R = utils.calculate_rotation_matrices(rotations)

# initialise
I    = np.random.random((classes, J, J))
#W    = np.empty((classes, rotations, pixels))
w    = np.ones((number_of_frames,))
b    = np.zeros((number_of_frames, 1))
logR = np.zeros((number_of_frames, classes, rotations))
P    = np.zeros((number_of_frames, classes, rotations))

# float32
print(f'writing reconstruction variables to: {recon_file}')
with h5py.File(recon_file, 'w') as f:
    f['rotation_matrices'] = R.astype(dtype)
    f['models'] = I.astype(dtype)
    f['model_voxel_size'] = dx_model
    f['model_xy_map'] = np.array(points_I).astype(dtype)
    f['fluence'] = w.astype(dtype)
    f['background_weights'] = b.astype(dtype)
    f['background'] = file_metadata['/entry_1/instrument_1/detector_1/background'][None, :].astype(dtype)
    f['logR'] = logR.astype(dtype)
    f['probability_matrix'] = P.astype(dtype)
    f['solid_angle_polarisation_factor'] = C.astype(dtype)
    f['iterations/iters'] = 0
    f.create_dataset('iterations/expectation_value', (config['iters'],), dtype=float, maxshape = (None,))
    f.create_dataset('iterations/log_likihood',      (config['iters'],), dtype=float, maxshape = (None,))
    f.create_dataset('iterations/beta',              (config['iters'],), dtype=float, maxshape = (None,))
    f.create_dataset('iterations/most_likely_class', (config['iters'], number_of_frames), dtype=float, maxshape = (None, number_of_frames))

# copy configuration file to output directory
shutil.copy2(sys.argv[1], out + '/config.py')
