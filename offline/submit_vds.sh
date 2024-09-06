#!/bin/bash

# call: sbatch submit_vds.sh <run_no>
# eg:   sbatch submit_vds.sh 2

source ../source_this_at_euxfel

ID=${EXP_ID}
PREFIX=${EXP_PREFIX}

#SBATCH --array=
#SBATCH --time=01:00:00
#SBATCH --export=ALL
#SBATCH -J vds
#SBATCH -o .vds-%4a-%j.out
#SBATCH -e .vds-%4a-%j.out
#SBATCH --partition=upex-beamtime
#SBATCH --reservation=upex_00${EXP_ID}
#####SBATCH --partition=upex

run=$(printf %.4d "${1}")
extra-data-make-virtual-cxi ${PREFIX}/raw/r${run} -o ${PREFIX}/scratch/vds/r${run}.cxi --exc-suspect-trains
