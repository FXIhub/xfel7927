#!/bin/bash

# call: ./submit_event_metadata.sh [get_event_metadata.py arguments]
# eg:   ./submit_event_metadata.sh
#       ./submit_event_metadata.sh -i <merged.cxi> -o <metadata.h5> -s Ery

source /etc/profile.d/modules.sh

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PARENT_DIR=$(dirname $SCRIPT_DIR)
source $PARENT_DIR/source_this_at_euxfel

cd $SCRIPT_DIR

sbatch <<EOT
#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH --export=ALL
#SBATCH -J metadata-${EXP_ID}
#SBATCH -o ${EXP_PREFIX}/scratch/log/metadata-${EXP_ID}-%j.out
#SBATCH -e ${EXP_PREFIX}/scratch/log/metadata-${EXP_ID}-%j.out
#SBATCH --partition=upex

# exit on first error
set -e

source /etc/profile.d/modules.sh
source $PARENT_DIR/source_this_at_euxfel

python get_event_metadata.py $@
EOT
