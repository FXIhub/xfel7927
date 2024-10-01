#!/bin/bash

# call: ./submit_is_hit.sh <run_no>
# eg:   ./submit_is_hit.sh 2

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
#SBATCH -J is_hit-${EXP_ID}
#SBATCH -o ${EXP_PREFIX}/scratch/log/is_hit-${EXP_ID}-%A-%a.out
#SBATCH -e ${EXP_PREFIX}/scratch/log/is_hit-${EXP_ID}-%A-%a.out
##SBATCH --partition=upex-beamtime
##SBATCH --reservation=upex_${EXP_ID}
#SBATCH --partition=upex

source /etc/profile.d/modules.sh
source $PARENT_DIR/source_this_at_euxfel

set -e

run=\${SLURM_ARRAY_TASK_ID}

# add is_hit
python add_is_hit.py \${run} -t 5

echo is_hit done
EOT
