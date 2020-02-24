#!/bin/bash
#script to gather data from scattered config files
#$ -l h_rt=08:00:00
#$ -cwd
#$ -l arch=intel*
#$ -l mem=6G
#$ -N gathloop_bv20G


#source /etc/profile.d/modules.sh add mpi/intel/openmpi/1.8.3
#module add mpi/intel/openmpi/1.8.3



echo "Run gather program on"
echo $HOSTNAME


#n=$1
n=702
i=669  
g=1 # step


    while [ $i -le $n  ]
    do

	./distribution -s=$i /fastdata/cs1mkg/smaug/sac3d_kinktest/3D_tube_vertical_np010404.out /fastdata/cs1mkg/smaug/sac3d_kinktest/3D_tube_vertical.out
	mv /fastdata/cs1mkg/smaug/sac3d_kinktest/3D_tube_vertical.out /fastdata/cs1mkg/smaug/sac3d_kinktest/3D_tube_vertical_$i.out
        i=`expr $i + $g`
        echo $i
        echo "/fastdata/cs1mkg/smaug/sac3d_kinktest/3D_tube_vertical.out_$i.out"


    done
