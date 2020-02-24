#!/bin/bash
#$ -l rmem=10G

#$ -pe mpi 8
##$ -l mem=32G

#$ -N S2D

echo "SAC will run on the following nodes"
echo $HOSTNAME
#cat $PE_HOSTFILE

#export LD_LIBRARY_PATH="/usr/local/packages/mpi/openmpi/2.0.1/intel-17.0.0/lib/:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/compilers_and_libraries_2017.0.098/linux/compiler/lib/intel64_lin:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/compilers_and_libraries_2017.0.098/linux/mpi/intel64/lib:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/compilers_and_libraries_2017.0.098/linux/mpi/mic/lib:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/compilers_and_libraries_2017.0.098/linux/ipp/lib/intel64:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/compilers_and_libraries_2017.0.098/linux/compiler/lib/intel64:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/compilers_and_libraries_2017.0.098/linux/mkl/lib/intel64:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/compilers_and_libraries_2017.0.098/linux/tbb/lib/intel64/gcc4.7:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/debugger_2017/iga/lib:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/debugger_2017/libipt/intel64/lib:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/compilers_and_libraries_2017.0.098/linux/daal/lib/intel64_lin:"$LD_LIBRARY_PATH
#export PATH="/usr/local/packages/mpi/openmpi/2.0.1/intel-17.0.0/bin/:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/compilers_and_libraries_2017.0.098/linux/bin/intel64:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/compilers_and_libraries_2017.0.098/linux/mpi/intel64/bin:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/debugger_2017/gdb/intel64_mic/bin:"$PATH





#source /etc/profile.d/modules.sh load mpi/openmpi/2.0.1/intel-17.0.0



echo $PATH
echo $LD_LIBRARY_PATH


module load mpi/openmpi/2.0.1/intel-17.0.0
 #module load mpi/openmpi/1.10.4/gcc-4.9.4-TESTING
#module load libs/CUDA/7.5.18/bina

#module load dev/intel-compilers/17.0.0
#export TIMECOUNTER=0
#ulimit -a

#which mpif90

#which mpirun





echo Starting vac now.
#source timeused
mpirun vac_sharc_mpi < vac.par
#mpirun -np 8 --hostfile sharchostfile8 vac_sharc_mpi < vac.par
#./vac < par/vac_OT.par
#source timeused


