#$ -cwd
#$ -V -S /bin/bash
#$ -N  conv_nn_TEMPFN_TEMPFW_TEMPH
#$ -pe smp 10
#$ -e err-TEMPFN_TEMPFW_TEMPH.log
#$ -o std-TEMPFN_TEMPFW_TEMPH.log
export PATH=~/miniconda3/bin:$PATH
hostname
mpirun -np 10 ./nd_mpi_TEMPFN_TEMPFW_TEMPH.py



