#!/bin/bash

# call: ./submit_events.sh <run_no>
# eg:   ./submit_events.sh 2

source /etc/profile.d/modules.sh

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PARENT_DIR=$(dirname $SCRIPT_DIR)
source $PARENT_DIR/source_this_at_euxfel

cd $SCRIPT_DIR

set -e

run=${1}

python make_events_file.py ${run} --hit_finding_mask hit_finding_mask.h5

# add pulse energy
python add_pulsedata.py ${run}

# add is_hit
python add_is_hit.py ${run} -t 4 --per_train

# calculate powder patterns for hits and non-hits
python save_powder_hits_nonhits.py ${run} 

echo events done

