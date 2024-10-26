#!/bin/bash

# call: ./submit_radial_profiles.sh <run_no> 

source /etc/profile.d/modules.sh

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PARENT_DIR=$(dirname $SCRIPT_DIR)
source $PARENT_DIR/source_this_at_euxfel

cd $SCRIPT_DIR

sbatch <<EOT
#!/bin/bash

#SBATCH --array=${1}
#SBATCH --time=03:00:00
#SBATCH --export=ALL
#SBATCH -J rad-${EXP_ID}
#SBATCH -o ${EXP_PREFIX}/scratch/log/rad-${EXP_ID}-%A-%a.out
#SBATCH -e ${EXP_PREFIX}/scratch/log/rad-${EXP_ID}-%A-%a.out
##SBATCH --partition=upex-beamtime
##SBATCH --reservation=upex_${EXP_ID}
#SBATCH --partition=upex

# exit on first error
set -e

source /etc/profile.d/modules.sh
source $PARENT_DIR/source_this_at_euxfel

run=\${SLURM_ARRAY_TASK_ID}
echo ${1} run = \${run} 

python save_radial_profiles.py \${run} -m r0551_mask.h5

echo rad done
EOT
