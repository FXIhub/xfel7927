#!/bin/bash

# call: sbatch submit_vds.sh <run_no>
# eg:   sbatch submit_vds.sh 2

source /etc/profile.d/modules.sh
source ../source_this_at_euxfel

sbatch <<EOT
#SBATCH --array=${1}
#SBATCH --time=01:00:00
#SBATCH --export=ALL
#SBATCH -J vds-${EXP_ID}-\${a}
#SBATCH -o ${EXP_PREFIX}/scratch/log/vds-${EXP_ID}-%A-%a.out
#SBATCH -e ${EXP_PREFIX}/scratch/log/vds-${EXP_ID}-%A-%a.out
###SBATCH --partition=upex-beamtime
###SBATCH --reservation=upex_${EXP_ID}
#SBATCH --partition=upex

source /etc/profile.d/modules.sh
source ../source_this_at_euxfel

ID=${EXP_ID}
PREFIX=${EXP_PREFIX}

run=$(printf %.4d "\${SLURM_ARRAY_TASK_ID}")
extra-data-make-virtual-cxi ${PREFIX}/raw/r\${run} -o ${PREFIX}/scratch/vds/r\${run}.cxi --exc-suspect-trains

echo vds done
EOT
