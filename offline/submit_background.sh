#!/bin/bash

# call: ./submit_background.sh <run_no>
# eg:   ./submit_background.sh 2

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
#SBATCH -J background-${EXP_ID}
#SBATCH -o ${EXP_PREFIX}/scratch/log/background-${EXP_ID}-%A-%a.out
#SBATCH -e ${EXP_PREFIX}/scratch/log/background-${EXP_ID}-%A-%a.out
#SBATCH --partition=upex-beamtime
#SBATCH --reservation=upex_${EXP_ID}
##SBATCH --partition=upex

source /etc/profile.d/modules.sh
source $PARENT_DIR/source_this_at_euxfel

set -e

run=\${SLURM_ARRAY_TASK_ID}

python add_background.py \${run} -n 32 

echo background done
EOT
