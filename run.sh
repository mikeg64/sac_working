#!/bin/bash
##$ -l h_rt=99:00:00
#$ -cwd
##$ -l arch=intel-x5650
##$ -l arch=intel-e5-2650v2
#$ -l mem=6G

#$ -N S2D



#module add compilers/intel/15.0.3

#export TIMECOUNTER=0

echo Starting vac now.
#source timeused
./vac < vac.par
#source timeused


