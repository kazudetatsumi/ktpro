#$ -cwd
#$ -V -S /bin/bash
#$ -N  conv_nn_test
#$ -pe smp 10
#$ -e err-job2.log
#$ -o std-job2.log
#$ -q all.q@nova3
export PATH=~/miniconda3/bin:$PATH
mpirun -np 10 ./nd_mpi_10_1_100.py



