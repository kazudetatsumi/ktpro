#$ -cwd
#$ -V -S /bin/bash
#$ -N phonon-num
#$ -pe mpi* 16
#$ -e err.log
#$ -o std.log

mpirun vasp541mpi


