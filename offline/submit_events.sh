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
#SBATCH --time=02:00:00
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

python make_events_file.py \${run} --hit_finding_mask hit_finding_mask.h5

# add pulse energy
python add_pulsedata.py \${run}

# add is_hit
python add_is_hit.py \${run} -t 4 --per_train

# calculate powder patterns for hits and non-hits
python powder.py \${run} -s is_hit True
python powder.py \${run} -s is_hit False

# calculate per-pixel powders 
# looks for files in scratch/powder/r<run>_powder_*.h5
python per_pixel_powder.py \${run} 

# estimate background strength per event
#python estimate_background.py \${run}

echo events done
EOT
