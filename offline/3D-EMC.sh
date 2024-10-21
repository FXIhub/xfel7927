#!/bin/bash
# Launch python code
start=`date +%s`

cd ${EXP_PREFIX}/scratch/3D-EMC/Ery

# remove cache
#rm -r cachedir/*
#rm *.pickle

#ln -s ${EXP_PREFIX}/scratch/saved_hits/Ery_classified_hits.cxi data.cxi

#qmin=0
#qmax=100000000 
#qmin_prob=30000000
#qmax_prob=90000000 
#mpx=128

# exit on first error
set -e

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PARENT_DIR=$(dirname $SCRIPT_DIR)

source /etc/profile.d/modules.sh
source $PARENT_DIR/source_this_at_euxfel
#source $PARENT_DIR/source_this_at_home

conda activate /home/amorgan/.conda/envs/EMC

#init_emc.py --qmin ${qmin} --qmax ${qmax} --mpx ${mpx}
#
#for i in {1..10}; do
#    logR.py --dc 8192 --ic 8192 --rc 8192 -M 20 --qmin ${qmin_prob} --qmax ${qmax_prob} 
#    calculate_probabilities.py --beta 0.01
#    update_I.py --ic 8192 --inversion_symmetry --qmax ${qmax} --p_thresh 0.00001 --mpx ${mpx} 
#done
#
#for i in {1..10}; do
#    logR.py --dc 8192 --ic 8192 --rc 8192 -M 20 --qmin ${qmin_prob} --qmax ${qmax_prob} 
#    calculate_probabilities.py --beta 0.05
#    update_I.py --ic 8192 --inversion_symmetry --qmax ${qmax} --p_thresh 0.00001 --mpx ${mpx} 
#done
#
qmin=0
qmax=300000000 
qmin_prob=30000000
qmax_prob=290000000 
mpx=256
#
#update_I.py --ic 8192 --inversion_symmetry --qmax ${qmax} --p_thresh 0.00001 --mpx ${mpx} 
#
#for i in {1..20}; do
#    logR.py --dc 8192 --ic 8192 --rc 8192 -M 20 --qmin ${qmin_prob} --qmax ${qmax_prob} 
#    calculate_probabilities.py --beta 0.05
#    update_I.py --ic 8192 --inversion_symmetry --qmax ${qmax} --p_thresh 0.00001 --mpx ${mpx} 
#done

#for i in {1..10}; do
#    logR.py --dc 8192 --ic 8192 --rc 8192 -M 20 --qmin ${qmin_prob} --qmax ${qmax_prob} 
#    calculate_probabilities.py --beta 0.1
#    update_I.py --ic 8192 --inversion_symmetry --qmax ${qmax} --p_thresh 0.00001 --mpx ${mpx} 
#done

#for i in {1..10}; do
#    logR.py --dc 8192 --ic 8192 --rc 8192 -M 20 --qmin ${qmin_prob} --qmax ${qmax_prob} 
#    calculate_probabilities.py --beta 0.2
#    update_I.py --ic 8192 --inversion_symmetry --qmax ${qmax} --p_thresh 0.00001 --mpx ${mpx} 
#done
#
#for i in {1..10}; do
#    logR.py --dc 8192 --ic 8192 --rc 8192 -M 20 --qmin ${qmin_prob} --qmax ${qmax_prob} 
#    calculate_probabilities.py --beta 0.4
#    update_I.py --ic 8192 --inversion_symmetry --qmax ${qmax} --p_thresh 0.00001 --mpx ${mpx} 
#done
#
#for i in {1..10}; do
#    logR.py --dc 8192 --ic 8192 --rc 8192 -M 20 --qmin ${qmin_prob} --qmax ${qmax_prob} 
#    calculate_probabilities.py --beta 1.0
#    update_I.py --ic 8192 --inversion_symmetry --qmax ${qmax} --p_thresh 0.00001 --mpx ${mpx} 
#done

update_I.py --ic 8192 --inversion_symmetry --qmax ${qmax} --p_thresh 0.01 --mpx ${mpx} 
update_I.py --ic 8192 --inversion_symmetry --merge_I --qmax ${qmax} --p_thresh 0.00001 --mpx ${mpx} 

end=`date +%s`
runtime=$((end-start))
echo 
echo runtime: "$runtime" seconds

