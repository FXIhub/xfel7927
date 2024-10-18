#!/bin/bash
# Launch python code
start=`date +%s`

# remove cache
#rm -r cachedir/*
#rm *.pickle

cd ${EXP_PREFIX}/scratch/3D-EMC/Ery

ln -s ${EXP_PREFIX}/scratch/saved_hits/Ery_classified_hits.cxi data.cxi

qmin=0
qmax=100000000 
qmin_prob=30000000
qmax_prob=90000000 
mpx=128

# exit on first error
set -e

init_emc.py --qmin ${qmin} --qmax ${qmax} --mpx ${mpx}

for i in {1..10}; do
    logR.py --dc 8192 --ic 8192 --rc 8192 -M 20 --qmin ${qmin_prob} --qmax ${qmax_prob} 
    calculate_probabilities.py --beta 0.05
    update_I.py --ic 8192 --inversion_symmetry --qmax ${qmax} --p_thresh 0.00001 --mpx ${mpx} 
done


end=`date +%s`
runtime=$((end-start))
echo 
echo runtime: "$runtime" seconds

