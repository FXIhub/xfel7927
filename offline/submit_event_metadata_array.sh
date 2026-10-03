#!/bin/bash

# one get_event_metadata.py job per run (slurm array), then a job that merges
# the per-run files once every run has succeeded
#
# call: ./submit_event_metadata_array.sh [sample] [merged cxi file] [output file]
# eg:   ./submit_event_metadata_array.sh
#       ./submit_event_metadata_array.sh Ery
#
# per-run files go to <output without .h5>_runs/r<run>.h5. If a run fails the
# merge job is cancelled; rerun that run and merge by hand:
#       python get_event_metadata.py --runs <run> -o <dir>/r<run>.h5
#       python get_event_metadata.py --merge <dir>/r*.h5 -o <output>

source /etc/profile.d/modules.sh

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PARENT_DIR=$(dirname $SCRIPT_DIR)
source $PARENT_DIR/source_this_at_euxfel

cd $SCRIPT_DIR

SAMPLE=${1:-Ery}
INPUT=${2:-${EXP_PREFIX}/scratch/saved_hits/${SAMPLE}_all_hits_no_mask.cxi}
OUTPUT=${3:-${EXP_PREFIX}/scratch/saved_hits/${SAMPLE}_event_metadata.h5}
RUN_DIR=${OUTPUT%.h5}_runs
mkdir -p $RUN_DIR

RUNS=$(python -c "from get_event_metadata import get_runs; print(','.join(str(r) for r in get_runs('${SAMPLE}')))")
echo "sample ${SAMPLE}: runs ${RUNS}"

ARRAY_ID=$(sbatch --parsable <<EOT
#!/bin/bash

#SBATCH --array=${RUNS}%64
#SBATCH --time=01:00:00
#SBATCH --export=ALL
#SBATCH -J metadata-${EXP_ID}
#SBATCH -o ${EXP_PREFIX}/scratch/log/metadata-${EXP_ID}-%A-%a.out
#SBATCH -e ${EXP_PREFIX}/scratch/log/metadata-${EXP_ID}-%A-%a.out
#SBATCH --partition=upex

# exit on first error
set -e

source /etc/profile.d/modules.sh
source $PARENT_DIR/source_this_at_euxfel

run=\${SLURM_ARRAY_TASK_ID}
python get_event_metadata.py --runs \${run} -i ${INPUT} -o ${RUN_DIR}/r\$(printf %04d \${run}).h5
EOT
)
echo "submitted array job ${ARRAY_ID}"

MERGE_ID=$(sbatch --parsable --dependency=afterok:${ARRAY_ID} --kill-on-invalid-dep=yes <<EOT
#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH --export=ALL
#SBATCH -J metadata-merge-${EXP_ID}
#SBATCH -o ${EXP_PREFIX}/scratch/log/metadata-merge-${EXP_ID}-%j.out
#SBATCH -e ${EXP_PREFIX}/scratch/log/metadata-merge-${EXP_ID}-%j.out
#SBATCH --partition=upex

set -e

source /etc/profile.d/modules.sh
source $PARENT_DIR/source_this_at_euxfel

python get_event_metadata.py --merge ${RUN_DIR}/r*.h5 -o ${OUTPUT}
EOT
)
echo "submitted merge job ${MERGE_ID} (runs after all of ${ARRAY_ID} succeed)"
