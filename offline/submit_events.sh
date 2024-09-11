#!/bin/bash

# call: sbatch submit_events.sh <run_no>
# eg:   sbatch submit_events.sh 2

source /etc/profile.d/modules.sh
source ../source_this_at_euxfel

sbatch <<EOT
#!/bin/bash

#SBATCH --array=${1}
#SBATCH --time=01:00:00
#SBATCH --export=ALL
#SBATCH -J events-${EXP_ID}-${a}
#SBATCH -o .events-${EXP_ID}-%4a-%j.out
#SBATCH -e .events-${EXP_ID}-%4a-%j.out
###SBATCH --partition=upex-beamtime
###SBATCH --reservation=upex_${EXP_ID}
#SBATCH --partition=upex

source /etc/profile.d/modules.sh
source ../source_this_at_euxfel

run=${SLURM_ARRAY_TASK_ID}
python make_events_file.py ${run} -n 32

# add pulse energy
python add_pulsedata.py ${run}

# add is_hit
python add_is_hit.py ${run}
EOT
