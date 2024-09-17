3# Offline data access

To access the offline computer cluster, login to the display gateway of the Maxwell cluster:

```
ssh -X <upex account>@max-display.desy.de
```

You will need X11 forwarding (using e.g. XQuartz) to be able to display graphics from the login node on your screen.

Don't run big jobs on the login nodes. The Maxwell cluster uses Slurm to submit analysis jobs to the cluster, use the `upex` partition for XFEL analysis. To submit a job and view its status:

```
sbatch -p upex <job command>
squeue -u <upex account>
```

To allocate an analysis node that you can use interactively, type:

```
srun -p upex -t 10:00:00 --pty $SHELL -i
```

You can find more useful Slurm commands here:

http://xray.bmc.uu.se/~filipe/admin/davinci_user.html

You can find more info about offline data analysis at EuXFEL here:

http://www.desy.de/~barty/cheetah/Cheetah/EuXFEL_data_analysis.html

## Clone this repository

To get access to the analysis scripts in this repository, create a folder and clone the repository:

```
git clone https://github.com/FXIhub/xfel7927.git

```

This clones the repository to the current folder using HTTPS, which is fine for read access. To get write access to the repository, setup your SSH keys properly and do instead:

```
git clone git@github.com:FXIhub/xfel7927.git
```

Before you start performing your analysis, don't forget to load the proper modules:

```
source xfel7927/source_this_at_euxfel
```

## Experiment folder

Once you login to the cluster you will start in your home directory. To go to the experiment directory, write:

```
cd /gpfs/exfel/exp/SPB/...
```

This contains the `proc` folder for processed runs, `raw` folder for raw data runs, `usr` for smaller user files like software and calibration files and `scratch` for larger data files. Please make a folder in `scratch` with your UPEX user name and place your analysis output for the experiment there:

```
cd scratch
mkdir <upex account>
```

For sharing files with other users, please make sure they have the correct permissions, which if you're being lazy means:

```
chmod 755 *
```

## .py vs .sh scripts
Most of python scripts (`*.py`) is this directory have a corresponding `submit_*.sh` script, that submits a script using SLURM which then runs the python script on the maxwell cluster and outputs the log files in the `/scratch/log` directory.

## VDS files
A useful way to look at the data is to use the virtual data set (VDS) feature of HDF5. These files can be generated with the `submit_vds.sh` script in this folder. Some runs should already be converted in the `/scratch/vds/` folder. These files have all the modules for a given train in the same dataset slice. Thus, one can use the simple HDF5/h5py API to access the data for a given frame.

```
$ h5ls -r r0018_vds_raw.h5
/                        Group
/INSTRUMENT              Group
/INSTRUMENT/SPB_DET_AGIPD1M-1 Group
/INSTRUMENT/SPB_DET_AGIPD1M-1/DET Group
/INSTRUMENT/SPB_DET_AGIPD1M-1/DET/image Group
/INSTRUMENT/SPB_DET_AGIPD1M-1/DET/image/data Dataset {16, 173008, 2, 512, 128}
/INSTRUMENT/SPB_DET_AGIPD1M-1/DET/image/trainId Dataset {173008}

$ h5ls -r r0018_vds_proc.h5
/                        Group
/INSTRUMENT              Group
/INSTRUMENT/SPB_DET_AGIPD1M-1 Group
/INSTRUMENT/SPB_DET_AGIPD1M-1/DET Group
/INSTRUMENT/SPB_DET_AGIPD1M-1/DET/image Group
/INSTRUMENT/SPB_DET_AGIPD1M-1/DET/image/data Dataset {16, 173008, 512, 128}
/INSTRUMENT/SPB_DET_AGIPD1M-1/DET/image/trainId Dataset {173008}
```
The extra dimension in the raw data has the gain (digital) data for the frame.

## Hitfinding
The virtual data sets can be used for hitinding with the `submit_events.sh` or `make_events_file.py` scripts in this folder. The `submit_events.sh` script produces an HDF5 file with lit pixel values for all shots, that should be placed in the `/scratch/events/` folder with read permissions to all users `chmod +r /scratch/events/r*_events.h5`. They can then be read by the `add_is_hit.py` script to determine a hit threshold and adds the field `is_hit` back into the events file. `add_is_hit.py` is automattically called by the `submit_events.sh` script. This script also sums all images in the run to produce whats called powder pattern, these can be found in `scratch/powder/r*_powder.h5`. See for expample, the events file in p7076:

```
usr/Shared/amorgan/xfel7927/offline$ ./submit_events.sh 50
Submitted batch job 9590883

scratch/events% h5ls -r events_r0050.h5
/                        Group
/cellId                  Dataset {997893}
/hitscore                Dataset {997893}
/is_hit                  Dataset {997893}
/is_miss                 Dataset {997893}
/litpixels               Dataset {997893}
/photons                 Dataset {997893}
/pulseId                 Dataset {997893}
/trainId                 Dataset {997893}

scratch/powder$ h5ls -r r0050_powder.h5 
/                        Group
/powder                  Dataset {16, 512, 128}
```
## Saving Hits
Once the the field `is_hit` has been added to the events file we can then save the hits in separate cxi file for faster read access and file transfer. This cxi file is similar to the vds files above, but with more information about the beam and the detector. To save hits at `scratch/saved_hits/r*_hits.cxi` run the `submit_cxi.sh` script with the run number and sample name as an arguement, e.g.:
```
usr/Shared/amorgan/xfel7927/offline$ ./submit_cxi.sh 50 Cubic DNA origami

scratch/saved_hits$ h5ls -r r0050_hits.cxi 
/                        Group
/entry_1                 Group
/entry_1/cellId          Dataset {3716}
/entry_1/data_1          Group
/entry_1/data_1/data     Soft Link {/entry_1/instrument_1/detector_1/data}
/entry_1/experiment_identifier Dataset {3716}
/entry_1/instrument_1    Group
/entry_1/instrument_1/data_1 Group
/entry_1/instrument_1/detector_1 Group
/entry_1/instrument_1/detector_1/background Dataset {16, 512, 128}
/entry_1/instrument_1/detector_1/data Dataset {3716, 16, 512, 128}
/entry_1/instrument_1/detector_1/lit_pixels Dataset {3716}
/entry_1/instrument_1/detector_1/mask Dataset {3716, 16, 512, 128}
/entry_1/instrument_1/detector_1/photon_counts Dataset {3716}
/entry_1/instrument_1/detector_1/pixel_area Dataset {SCALAR}
/entry_1/instrument_1/detector_1/powder Dataset {16, 512, 128}
/entry_1/instrument_1/detector_1/x_pixel_size Dataset {SCALAR}
/entry_1/instrument_1/detector_1/xyz_map Dataset {3, 16, 512, 128}
/entry_1/instrument_1/detector_1/y_pixel_size Dataset {SCALAR}
/entry_1/instrument_1/name Dataset {SCALAR}
/entry_1/instrument_1/source_1 Group
/entry_1/instrument_1/source_1/photon_energy Dataset {3716}
/entry_1/instrument_1/source_1/photon_wavelength Dataset {3716}
/entry_1/instrument_1/source_1/pulse_energy Dataset {3716}
/entry_1/name            Dataset {SCALAR}
/entry_1/pulseId         Dataset {3716}
/entry_1/sample_1        Group
/entry_1/sample_1/name   Dataset {SCALAR}
/entry_1/trainId         Dataset {3716}
```

After `submit_cxi.sh` has called the `make_cxi_file.py` script and written the cxi file, it then calls the `cxi_to_emc.py` script, which writes the photon counts in a format that Dragonfly can read. These files are located in `scratch/emc/r*.emc`.
