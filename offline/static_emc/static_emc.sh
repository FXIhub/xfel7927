#!/bin/bash

set -e

conda activate /home/amorgan/.conda/envs/EMC

# make config file
run=$(printf %.4d "${1}")
static_dir=${EXP_PREFIX}scratch/static_emc
cxi_dir=${EXP_PREFIX}scratch/saved_hits
geom=${EXP_PREFIX}usr/geom/motor_p7076_from_4462.geom
mask=${EXP_PREFIX}scratch/det/static_emc_mask.h5

config=$(python make_config.py ${run} --cxi_dir ${cxi_dir} --output_dir ${static_dir} --mask ${mask} --geom ${geom})

# run init
python static_emc_init.py ${config}

# run emc 
mpirun -np 4 python static_emc.py ${config}

# plot results 
python plot_iters.py ${config}

