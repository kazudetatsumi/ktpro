#$ -cwd
#$ -V -S /bin/bash
#$ -N  two-layer-nn
#$ -pe x36 36
#$ -e err-NM.log
#$ -o std-NM.log
export PATH=~/miniconda3/bin:$PATH
mpirun -np 36 ./test2_mpi.py



