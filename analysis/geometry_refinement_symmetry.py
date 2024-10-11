import argparse

parser = argparse.ArgumentParser(description = 
r"""use a powder of bright hits to refine geometry using extra-geom

outputs a crystfel geometry file in the local directory
""")
parser.add_argument('powder_file', type=str, help='file name of powder file (*.h5) e.g. /gpfs/exfel/exp/SQS/202302/p003004/scratch/amorgan/powder_for_geometry_refinement_r456-r476.h5')
parser.add_argument('mask_file', type=str, help='file name of mask file (*.h5) (must contain entry_1/good_pixels or /good_pixels) e.g. /gpfs/exfel/exp/SQS/202302/p003004/scratch/det/backgroundpixel_mask_r0450.h5')
parser.add_argument('geom_file', type=str, help='file name of crystfel geometry file (*.geom)')
parser.add_argument('-q', '--refine_quads', action='store_true', help='refine quadrant positions in addition to beam centre')
parser.add_argument('-z', '--detector_distance', type=float, default=70e-2, help='detector distance will be written into geometry file')
parser.add_argument('-i', '--index', type=int, help='load a specific index in the powder file')
parser.add_argument('-d', '--data_path', type=str, default = 'data', help='the h5 path of the dataset in powder file')
parser.add_argument('-o', '--write_output', action = 'store_true', help='write output to file, otherwise the default behaviour is to show a plot')
args = parser.parse_args()
    
args.output_crystfel = f'crystfel_geometry_file.geom'


import numpy as np
import h5py
import pickle
import os.path as op


def rad_av(ar, mask, xyz):
    """
    xyz should be scaled of integer binning
    """
    #xyz = make_xyz(corner, basis, ar.shape)
    r = np.rint((xyz[0]**2 + xyz[1]**2)**0.5).astype(int)
    
    rsum = np.bincount(r.ravel(), (mask*ar).ravel())
    
    rcounts = np.bincount(r.ravel(), mask.ravel())
    
    rcounts[rcounts==0] = 1
    
    rav = rsum / rcounts
    
    # broacast
    ar_rav = np.zeros_like(ar)
    
    ar_rav[:] = rav[r.ravel()].reshape(ar.shape)
    
    return rav, ar_rav


import extra_geom


# load powder
with h5py.File(args.powder_file, 'r') as f:
    if args.index != None :
        powder = f[args.data_path][args.index]
    else :
        powder = f[args.data_path][()]

powder = powder.astype(float)

# load mask
with h5py.File(args.mask_file, 'r') as f:
    if 'entry_1' in f :
        mask = f['entry_1/good_pixels'][()]
    else :
        mask = f['good_pixels'][()]

geom = extra_geom.AGIPD_1MGeometry.from_crystfel_geom(args.geom_file)

# make x, y, z the first axis
xyz = np.transpose(geom.get_pixel_positions(), (3, 0, 1, 2))

r = (xyz[0]**2 + xyz[1]**2)**.5
rmask = (r > 0.0) * (r < 0.02)
rmask *= mask

# scale xyz for integer binning (approximately convert to pixel units)
xyz /= 236e-6



r_powder, powder_rav = rad_av(powder, mask, xyz)
im     = geom.position_modules(mask * powder)[0]
im_rav = geom.position_modules(powder_rav)[0]

im0 = im.copy()
im_rav0 = im_rav.copy()
powder_rav0 = powder_rav.copy()


def fill_image(im_in):
    t = [min(i, j) for i, j in zip(im.shape, im_in.shape)]
    im.fill(0)
    im[:t[0], :t[1]] = im_in[:t[0], :t[1]]
    return im

# first refine beam centre with 300um step size
ermin      = np.inf
step_size  = 400e-6
small_step_size  = 100e-6
shift_min  = np.zeros((16, 2), dtype=float)
shift      = np.zeros((16, 2), dtype=float)
ims = []

# if the user only wants beam centre refinement then repeat with smaller step size
if args.refine_quads == False :
    step_sizes = [step_size, small_step_size]
else :
    step_sizes = [step_size,]

for step in step_sizes:
    for i in range(1000):
        done = True
        for dim in range(2):
            # take a step and see if the error improves
            for j in [-1, 1]:
                shift[:]       = shift_min
                shift[:, dim] += j * step
                
                # this seems dodgy
                geom2 = geom.offset(shift)
                #geom2 = get_geom(quad_pos_pix + shift)
                 
                xyz = np.transpose(geom2.get_pixel_positions(), (3, 0, 1, 2)) / 236e-6
                
                # calculate error
                r_powder, powder_rav = rad_av(powder, rmask, xyz)
                err = np.sum(rmask * (powder - powder_rav)**2)
                 
                # take step if error improves
                if err < ermin :
                    shift_min[:] = shift             
                    ermin     = err
                    print(i, ermin, shift_min[0,0])
                    done = False

        if done :
            break

if args.refine_quads :
    # now refine quadrant positions (but not Q2) with 100um step size
    ermin      = np.inf
    step       = small_step_size
    q_shift_min  = np.zeros((16, 2), dtype=float)
    quadrant_slices = [np.s_[:4], np.s_[4:8], np.s_[8:12], np.s_[12:16]]
    # add global offset
    q_shift_min += shift_min
    ims = []

    for i in range(1000):
        done = True
        for Q in [0, 1, 2, 3]:
            for dim in range(2):
                # take a step and see if the error improves
                for j in [-1, 1]:
                    shift          = q_shift_min.copy()
                    shift[quadrant_slices[Q], dim] += j * step
                    
                    # this seems dodgy
                    geom2 = geom.offset(shift)
                    #geom2 = get_geom(quad_pos_pix + shift)
                     
                    xyz = np.transpose(geom2.get_pixel_positions(), (3, 0, 1, 2)) / 236e-6
                    
                    # calculate error
                    r_powder, powder_rav = rad_av(powder, rmask, xyz)
                    err = np.sum(rmask * (powder - powder_rav)**2)
                    
                    #ims.append(geom2.position_modules(rmask * (powder - powder_rav))[0])
                     
                    # take step if error improves
                    if err < ermin :
                        q_shift_min[:] = shift
                        ermin     = err
                        print(i, ermin, q_shift_min[::4])
                        ims.append(geom2.position_modules(powder)[0])
                        done = False
        if done :
            break
else :
    q_shift_min = shift_min


geom2 = geom.offset(q_shift_min)

xyz   = np.transpose(geom2.get_pixel_positions(), (3, 0, 1, 2)) / 236e-6
r_powder, powder_rav = rad_av(powder, rmask, xyz)


# add detector distance
geom2 = geom2.offset((0, 0, args.detector_distance))

# plot powder and initial powder azimuthally averaged 
if not args.write_output :
    import pyqtgraph as pg
    im_rav = geom.position_modules(powder_rav)[0]
    pg.show(np.array([im0, im_rav0, im0, im_rav])**0.2)

    if len(ims) > 0 :
        #ims.append(im)
        min_ss = min([a.shape[0] for a in ims])
        min_fs = min([a.shape[1] for a in ims])
        ims = np.array([a[:min_ss, :min_fs] for a in ims])
        pg.show(np.abs(np.array(ims)))

# has to be in mm due to bug I think
#quad_pos_pix = 1e3 * (quad_pos_pix + shift)

print('refined quad positions in um:')
print(q_shift_min[::4] * 1e6)

# save
import pickle
from pathlib import Path
fnam = f'quad_positions_{Path(args.powder_file).stem}.pickle'
pickle.dump( q_shift_min[::4], open(fnam, 'wb'))

if args.write_output :
    # write geom_<run_label>.py to scratch/det/
    # output crystfel geom file
    print('writing to', args.output_crystfel)
    geom2.write_crystfel_geom(args.output_crystfel)

print('done!')

# plot powder and ring
if not args.write_output :
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    from mpl_toolkits.axes_grid1 import make_axes_locatable


    #for d in [powder, np.abs(powder - powder_rav)] :
    for d in [powder, ] :
        ax = geom2.plot_data(mask * d, axis_units='m', vmin=0, vmax=d.max() / 10)
        #ax = geom2.plot_data_hexes(mask * d, vmin=0, vmax=d.max() / 10, colorbar=True)

        circ = Circle((0, 0), 0.0065, fill=False, color='r', linewidth=2, alpha=0.3)
        ax.add_patch(circ)

        circ = Circle((0, 0), 0.009, fill=False, color='r', linewidth=2, alpha=0.3)
        ax.add_patch(circ)

        ax.set_title(f'refined geometry from powder run')

        # make colorbar a good size
        im = ax.get_images()[0]
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        #plt.colorbar(im, cax=cax)



    fig2, ax = plt.subplots()

    xyz1 = np.transpose(geom.get_pixel_positions(), (3, 0, 1, 2)) / 236e-6
    xyz2 = np.transpose(geom2.get_pixel_positions(), (3, 0, 1, 2)) / 236e-6

    # calculate error
    r_powder1, powder_rav1 = rad_av(powder, rmask, xyz1)
    r_powder2, powder_rav2 = rad_av(powder, rmask, xyz2)

    r1 =  np.arange(r_powder1.shape[0]) * 236e-6
    r2 =  np.arange(r_powder2.shape[0]) * 236e-6
    ax.plot(r1, r_powder1, label='original radial average')
    ax.plot(r2, r_powder2, label='refined geometry radial average')
    ax.set_yscale('log')
    ax.set_xlim([0, 0.04])
    ax.legend()
    plt.show()

