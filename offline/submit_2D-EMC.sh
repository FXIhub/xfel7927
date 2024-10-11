#!/bin/bash

# need to run init_2D-EMC.sh first

# call: ./submit_2D-EMC.sh <recon directory>
# eg:   ./submit_2D-EMC.sh /gpfs/exfel/exp/SPB/202405/p007927/scratch/2D-EMC/Ery


# this is just to get environment variables defined in source_this_at_euxfel
source /etc/profile.d/modules.sh
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PARENT_DIR=$(dirname $SCRIPT_DIR)
source $PARENT_DIR/source_this_at_euxfel


sbatch <<EOT
#!/bin/bash

#SBATCH --time=010:00:00
#SBATCH --export=ALL
#SBATCH -J 2DEMC-${EXP_ID}
#SBATCH -o ${EXP_PREFIX}/scratch/log/2DEMC-${EXP_ID}-%A-%a.out
#SBATCH -e ${EXP_PREFIX}/scratch/log/2DEMC-${EXP_ID}-%A-%a.out
##SBATCH --partition=upex-beamtime
##SBATCH --reservation=upex_${EXP_ID}
#SBATCH --constraint="GPUx4&A100"
#SBATCH --partition=allgpu

# exit on first error
set -e

./run_2D-EMC.sh \${1}

echo 2DEMC done
EOT

