#!/bin/bash

# call: ./run_2D-EMC.sh <recon folder>
# eg:   ./run_2D-EMC.sh /gpfs/exfel/exp/SPB/202405/p007927/scratch/2D-EMC/Ery

# exit on first error
set -e

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PARENT_DIR=$(dirname $SCRIPT_DIR)

source /etc/profile.d/modules.sh
source $PARENT_DIR/source_this_at_euxfel
#source $PARENT_DIR/source_this_at_home

conda activate /home/amorgan/.conda/envs/EMC

cd $SCRIPT_DIR

cd 2D-EMC

mpirun -n 4 python 2D-EMC.py ${1}

#pdfunite ${1}/*_*.pdf ${1}/recon.pdf
