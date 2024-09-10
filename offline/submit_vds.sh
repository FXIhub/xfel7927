#!/bin/bash

# call: sbatch submit_vds.sh <run_no>
# eg:   sbatch submit_vds.sh 2

#SBATCH --array=${1}
#SBATCH --time=01:00:00
#SBATCH --export=ALL
#SBATCH -J vds-${EXP_ID}-%4a-%j
#SBATCH -o .vds-${EXP_ID}-%4a-%j.out
#SBATCH -e .vds-${EXP_ID}-%4a-%j.out
###SBATCH --partition=upex-beamtime
###SBATCH --reservation=upex_${EXP_ID}
#SBATCH --partition=upex

source /etc/profile.d/modules.sh
source ../source_this_at_euxfel

ID=${EXP_ID}
PREFIX=${EXP_PREFIX}


run=$(printf %.4d "${1}")
extra-data-make-virtual-cxi ${PREFIX}/raw/r${run} -o ${PREFIX}/scratch/vds/r${run}.cxi --exc-suspect-trains
