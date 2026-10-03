#!/bin/bash

# write every frame of each run to a publication cxi file (one job per run):
# per-event metadata from the run's VDS, then make_publication_run_cxi.py
#
# call: ./submit_publication_runs.sh <runs>
# eg:   ./submit_publication_runs.sh 642-644
#       ./submit_publication_runs.sh 99-101,331-333,478-480
#
# output: ${EXP_PREFIX}/scratch/publication/r<run>_metadata.h5
#         ${EXP_PREFIX}/scratch/publication/p007927_r<run>_gas_background.cxi

source /etc/profile.d/modules.sh

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PARENT_DIR=$(dirname $SCRIPT_DIR)
source $PARENT_DIR/source_this_at_euxfel

cd $SCRIPT_DIR
OUT=${EXP_PREFIX}/scratch/publication
mkdir -p $OUT

sbatch <<EOT
#!/bin/bash

#SBATCH --array=${1}
#SBATCH --time=04:00:00
#SBATCH --export=ALL
#SBATCH -J publication-${EXP_ID}
#SBATCH -o ${EXP_PREFIX}/scratch/log/publication-${EXP_ID}-%A-%a.out
#SBATCH -e ${EXP_PREFIX}/scratch/log/publication-${EXP_ID}-%A-%a.out
#SBATCH --partition=upex

# exit on first error
set -e

source /etc/profile.d/modules.sh
source $PARENT_DIR/source_this_at_euxfel

run=\${SLURM_ARRAY_TASK_ID}
r=\$(printf %04d \${run})

python get_event_metadata.py -i ${EXP_PREFIX}/scratch/vds/r\${r}.cxi --runs \${run} -o ${OUT}/r\${r}_metadata.h5
python make_publication_run_cxi.py \${run} -n 32

echo publication done
EOT
