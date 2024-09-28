#!/bin/bash

# call: ./submit_intensity.sh 

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
#SBATCH -J intensity-${EXP_ID}
#SBATCH -o ${EXP_PREFIX}/scratch/log/intensity-${EXP_ID}-%A-%a.out
#SBATCH -e ${EXP_PREFIX}/scratch/log/intensity-${EXP_ID}-%A-%a.out
#SBATCH --partition=upex-beamtime
#SBATCH --reservation=upex_${EXP_ID}
##SBATCH --partition=upex

# exit on first error
set -e

source /etc/profile.d/modules.sh
source $PARENT_DIR/source_this_at_euxfel

python show_peak_intensities.py 

echo intensity done
EOT
