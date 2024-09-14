#!/bin/bash

# call: ./submit_cxi.sh <run_no> -s <sample name>
# eg:   ./submit_events.sh 2 DNA

source /etc/profile.d/modules.sh
source ../source_this_at_euxfel

sbatch <<EOT
#!/bin/bash

#SBATCH --array=${1}
#SBATCH --time=01:00:00
#SBATCH --export=ALL
#SBATCH -J cxi-${EXP_ID}-\${a}
#SBATCH -o ${EXP_PREFIX}/scratch/log/cxi-${EXP_ID}-%A-%a.out
#SBATCH -e ${EXP_PREFIX}/scratch/log/cxi-${EXP_ID}-%A-%a.out
###SBATCH --partition=upex-beamtime
###SBATCH --reservation=upex_${EXP_ID}
#SBATCH --partition=upex

# exit on first error
set -e

source /etc/profile.d/modules.sh
source ../source_this_at_euxfel

run=\${SLURM_ARRAY_TASK_ID}
echo ${1} run = \${run}   sample = ${2}

python make_cxi_file.py \${run} -s ${2}

# add background
python add_background_cxi.py \${run}

# add powder of hits
python add_powder_cxi.py \${run}

echo cxi done

# write emc file
python cxi_to_emc.py \${run}

echo emc done
EOT
