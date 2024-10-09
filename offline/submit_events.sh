#!/bin/bash

# call: ./submit_events.sh <run_no>
# eg:   ./submit_events.sh 2

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
#SBATCH -J events-${EXP_ID}
#SBATCH -o ${EXP_PREFIX}/scratch/log/events-${EXP_ID}-%A-%a.out
#SBATCH -e ${EXP_PREFIX}/scratch/log/events-${EXP_ID}-%A-%a.out
##SBATCH --partition=upex-beamtime
##SBATCH --reservation=upex_${EXP_ID}
#SBATCH --partition=upex

source /etc/profile.d/modules.sh
source $PARENT_DIR/source_this_at_euxfel

set -e

run=\${SLURM_ARRAY_TASK_ID}

python make_events_file.py \${run} -n 32 --masks r0572_mask.h5 hit_finding_mask.h5

# add pulse energy
python add_pulsedata.py \${run}

# add is_hit
python add_is_hit.py \${run} -t 3 --per_train

echo events done
EOT
