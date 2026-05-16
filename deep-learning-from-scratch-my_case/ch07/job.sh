#$ -cwd
#$ -V -S /bin/bash
#$ -N  conv_nn_test2_ndK2
#$ -pe smp 1
#$ -e err-test2_ndK2.log
#$ -o std-test2_ndK2.log
export PATH=~/miniconda3/bin:$PATH
time mpirun -np 1 ./test2_ndK2.py



