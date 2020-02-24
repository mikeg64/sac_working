#!/bin/bash
#script to gather data from scattered config files
#$ -l h_rt=08:00:00
#$ -cwd
#$ -l arch=intel*
#$ -l mem=6G
#$ -N gather_bv20G


#source /etc/profile.d/modules.sh add mpi/intel/openmpi/1.8.3
#module add mpi/intel/openmpi/1.8.3



echo "Run gather program on"
echo $HOSTNAME


data/distribution /fastdata/cs1mkg/sac/p5b0_1_bv20g/3D_spic_128_bv20G_np080101.out /fastdata/cs1mkg/sac/p5b0_1_bv20g/3D_spic_128_bv20G.out
