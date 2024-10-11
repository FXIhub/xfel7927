#!/bin/bash

# call: ./init_2D-EMC.sh <config file> <recon folder>
# eg:   ./init_2D-EMC.sh config_Ery_sized.py Ery_sized

# the config files are in the 2D-EMC folder
# the recon folder is relative to {PREFIX}/scratch/2D-EMC

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PARENT_DIR=$(dirname $SCRIPT_DIR)

source /etc/profile.d/modules.sh
source $PARENT_DIR/source_this_at_euxfel
#source $PARENT_DIR/source_this_at_home

cd $SCRIPT_DIR

cd 2D-EMC

if [ $# -eq 2 ]; then 
    python init.py ${1} ${EXP_PREFIX}/scratch/2D-EMC/${2}
else 
    echo no command line arguments! $#
fi

