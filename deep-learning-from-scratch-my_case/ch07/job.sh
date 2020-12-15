#$ -cwd
#$ -V -S /bin/bash
#$ -N  conv_nn_test3_mpi
#$ -pe smp 10
#$ -e err-test3_nd_mpi.log
#$ -o std-test3_nd_mpi.log
export PATH=~/miniconda3/bin:$PATH
mpirun -np 10 ./test2_nd_mpi.py



