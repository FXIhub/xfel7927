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
#SBATCH -J cxi-${EXP_ID}
#SBATCH -o ${EXP_PREFIX}/scratch/log/cxi-${EXP_ID}-%A-%a.out
#SBATCH -e ${EXP_PREFIX}/scratch/log/cxi-${EXP_ID}-%A-%a.out
#SBATCH --partition=upex-beamtime
#SBATCH --reservation=upex_${EXP_ID}
##SBATCH --partition=upex

# exit on first error
set -e

source /etc/profile.d/modules.sh
source $PARENT_DIR/source_this_at_euxfel

run=\${SLURM_ARRAY_TASK_ID}
echo ${1} run = \${run}   sample = ${2}

python make_cxi_file.py \${run} -s ${2}

# add background
python add_background_cxi.py \${run}

echo cxi done

# write emc file
python cxi_to_emc.py \${run}

echo emc done
EOT
