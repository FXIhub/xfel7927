#!/bin/bash

# call: ./submit_cxi.sh <run_no> -s <sample name>
# eg:   ./submit_events.sh 2 DNA

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
#SBATCH -J static_emc-${EXP_ID}
#SBATCH -o ${EXP_PREFIX}/scratch/log/static_emc-${EXP_ID}-%A-%a.out
#SBATCH -e ${EXP_PREFIX}/scratch/log/static_emc-${EXP_ID}-%A-%a.out
###SBATCH --partition=upex-beamtime
###SBATCH --reservation=upex_${EXP_ID}
#SBATCH --partition=upex

# exit on first error
set -e

source /etc/profile.d/modules.sh
source $PARENT_DIR/source_this_at_euxfel

cd static_emc

run=\${SLURM_ARRAY_TASK_ID}

./static_emc/static_emc.sh \${run}

echo static_emc done
EOT

