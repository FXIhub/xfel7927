#!/bin/bash

# call: ./submit_vds.sh <run_no>
# eg:   ./submit_vds.sh 2

source /etc/profile.d/modules.sh

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PARENT_DIR=$(dirname $SCRIPT_DIR)
source $PARENT_DIR/source_this_at_euxfel

sbatch <<EOT
#!/bin/bash

#SBATCH --array=${1}
#SBATCH --time=01:00:00
#SBATCH --export=ALL
#SBATCH -J vds-${EXP_ID}-\${a}
#SBATCH -o ${EXP_PREFIX}/scratch/log/vds-${EXP_ID}-%A-%a.out
#SBATCH -e ${EXP_PREFIX}/scratch/log/vds-${EXP_ID}-%A-%a.out
###SBATCH --partition=upex-beamtime
###SBATCH --reservation=upex_${EXP_ID}
#SBATCH --partition=upex

set -e 

source /etc/profile.d/modules.sh
source $PARENT_DIR/source_this_at_euxfel

run=\$(printf %.4d "\${SLURM_ARRAY_TASK_ID}")
extra-data-make-virtual-cxi ${EXP_PREFIX}/proc/r\${run} -o ${EXP_PREFIX}/scratch/vds/r\${run}.cxi --exc-suspect-trains

echo vds done
EOT
