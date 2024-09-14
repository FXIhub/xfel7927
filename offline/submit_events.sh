#!/bin/bash

# call: ./submit_events.sh <run_no>
# eg:   ./submit_events.sh 2

source /etc/profile.d/modules.sh
source ../source_this_at_euxfel

sbatch <<EOT
#!/bin/bash

#SBATCH --array=${1}
#SBATCH --time=01:00:00
#SBATCH --export=ALL
#SBATCH -J events-${EXP_ID}-\${a}
#SBATCH -o ${EXP_PREFIX}/scratch/log/events-${EXP_ID}-%A-%a.out
#SBATCH -e ${EXP_PREFIX}/scratch/log/events-${EXP_ID}-%A-%a.out
###SBATCH --partition=upex-beamtime
###SBATCH --reservation=upex_${EXP_ID}
#SBATCH --partition=upex

source /etc/profile.d/modules.sh
source ../source_this_at_euxfel

set -e

run=\${SLURM_ARRAY_TASK_ID}
python make_events_file.py \${run} -n 32

# add pulse energy
python add_pulsedata.py \${run}

# add is_hit
python add_is_hit.py \${run}

echo events done
EOT
