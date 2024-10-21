#!/bin/bash
# Launch python code
start=`date +%s`

# exit on first error
set -e

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PARENT_DIR=$(dirname $SCRIPT_DIR)

source /etc/profile.d/modules.sh
source $PARENT_DIR/source_this_at_euxfel
#source $PARENT_DIR/source_this_at_home

conda activate /home/amorgan/.conda/envs/EMC

cd ${EXP_PREFIX}/scratch/3D-EMC/Ery

# conver emc output to phasing input
emc_to_phasing_input.py < merged_intensity.pickle --beamstop 32 --cut_edges --soft_edges --size 300 -o intensity_phasing.pickle

# phase intensity_phasing and produce output-phasing.pickle
phase.py -s -b --no-centre --inversion_symmetry --shrink 1.5 0.4 20 --iters 100DM 1000ERA 100DM 1000ERA 100DM 1000ERA < intensity_phasing.pickle > output-phasing.pickle

python save_compressed_phase.py

end=`date +%s`
runtime=$((end-start))
echo 
echo runtime: "$runtime" seconds


