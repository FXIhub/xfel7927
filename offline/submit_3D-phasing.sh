#!/bin/bash

# need to run submit_3D-EMC.sh first

# call: ./submit_3D-phasing.sh <recon directory>
# eg:   ./submit_3D-phasing.sh /gpfs/exfel/exp/SPB/202405/p007927/scratch/2D-EMC/Ery


# this is just to get environment variables defined in source_this_at_euxfel
source /etc/profile.d/modules.sh
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PARENT_DIR=$(dirname $SCRIPT_DIR)
source $PARENT_DIR/source_this_at_euxfel


sbatch <<EOT
#!/bin/bash

#SBATCH --time=010:00:00
#SBATCH --export=ALL
#SBATCH -J 3Dphasing-${EXP_ID}
#SBATCH -o ${EXP_PREFIX}/scratch/log/3Dphasing-${EXP_ID}.out
#SBATCH -e ${EXP_PREFIX}/scratch/log/3Dphasing-${EXP_ID}.out
##SBATCH --partition=upex-beamtime
##SBATCH --reservation=upex_${EXP_ID}
#SBATCH --constraint="A100"
#SBATCH --partition=allgpu

# exit on first error
set -e

./3D-phasing.sh

echo 3Dphasing done
EOT



