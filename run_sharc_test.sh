#!/bin/bash
##$ -l h_rt=99:00:00
##$ -cwd
##$ -l arch=intel-x5650
##$ -l arch=intel-e5-2650v2
##$ -l rmem=64G
##$ -pe openmpi-ib 8
#$ -pe mpi 16

#$ -N S2D

echo "SAC will run on the following nodes"
cat $PE_HOSTFILE > sharchostfile

#export LD_LIBRARY_PATH="/usr/local/packages/mpi/openmpi/2.0.1/intel-17.0.0/lib/:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/compilers_and_libraries_2017.0.098/linux/compiler/lib/intel64_lin:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/compilers_and_libraries_2017.0.098/linux/mpi/intel64/lib:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/compilers_and_libraries_2017.0.098/linux/mpi/mic/lib:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/compilers_and_libraries_2017.0.098/linux/ipp/lib/intel64:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/compilers_and_libraries_2017.0.098/linux/compiler/lib/intel64:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/compilers_and_libraries_2017.0.098/linux/mkl/lib/intel64:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/compilers_and_libraries_2017.0.098/linux/tbb/lib/intel64/gcc4.7:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/debugger_2017/iga/lib:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/debugger_2017/libipt/intel64/lib:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/compilers_and_libraries_2017.0.098/linux/daal/lib/intel64_lin:"$LD_LIBRARY_PATH
#export PATH="/usr/local/packages/mpi/openmpi/2.0.1/intel-17.0.0/bin/:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/compilers_and_libraries_2017.0.098/linux/bin/intel64:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/compilers_and_libraries_2017.0.098/linux/mpi/intel64/bin:/usr/local/packages/dev/intel-ps-xe-ce/2017.0/binary/debugger_2017/gdb/intel64_mic/bin:"$PATH





#source /etc/profile.d/modules.sh load mpi/openmpi/2.0.1/intel-17.0.0




module load mpi/openmpi/2.0.1/intel-17.0.0

#export TIMECOUNTER=0


which mpif90

which mpirun


echo $PATH
echo $LD_LIBRARY_PATH

cat sharchostfile



echo Starting vac now.
#source timeused
mpirun vac_sharc_mpi < vac.par
#mpirun -np 8 --hostfile sharchostfile8 vac_sharc_mpi < vac.par

#source timeused


