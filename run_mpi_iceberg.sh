#!/bin/bash
#$ -l h_rt=99:00:00
#$ -cwd
##$ -q cstest.q
##$ -l arch=intel-x5650
#$ -l arch=intel-e5-2650v2
##$ -l mem=58G
##$ -l rmem=58G
##$ -P cstest
#$ -pe openmpi-ib 16

#$ -N S2D

echo "SAC will run on the following nodes"
cat $PE_HOSTFILE > samplehostfile


module add mpi/intel/openmpi/1.8.3

export TIMECOUNTER=0

echo Starting vac now.
source timeused
mpirun vac_iceberg_mpi < vac.par
source timeused


