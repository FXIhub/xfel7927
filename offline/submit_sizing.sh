
#!/bin/bash

# call: ./submit_sizing.sh <run_no> 
# eg:   ./submit_sizing.sh 2 

source /etc/profile.d/modules.sh

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PARENT_DIR=$(dirname $SCRIPT_DIR)
source $PARENT_DIR/source_this_at_euxfel

cd $SCRIPT_DIR

sbatch <<EOT
#!/bin/bash

#SBATCH --array=${1}
#SBATCH --time=01:00:00
#SBATCH --export=ALL
#SBATCH -J sizing-${EXP_ID}
#SBATCH -o ${EXP_PREFIX}/scratch/log/sizing-${EXP_ID}-%A-%a.out
#SBATCH -e ${EXP_PREFIX}/scratch/log/sizing-${EXP_ID}-%A-%a.out
##SBATCH --partition=upex-beamtime
##SBATCH --reservation=upex_${EXP_ID}
#SBATCH --partition=upex

# exit on first error
set -e

source /etc/profile.d/modules.sh
source $PARENT_DIR/source_this_at_euxfel

# this is to get mpi4py (overkill)
conda activate /home/amorgan/.conda/envs/EMC

run=\${SLURM_ARRAY_TASK_ID}
echo ${1} run = \${run} 

#mpirun -n 32 python add_sizing_cxi.py \${run} --N_sizes 32 --N_angles 32
python add_sizing_cxi.py \${run} --N_sizes 64 --N_angles 64

echo sizing done
EOT
