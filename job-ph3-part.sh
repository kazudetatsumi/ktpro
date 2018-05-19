#$ -cwd
#$ -V -S /bin/bash
#$ -N ph3gp_noig
#$ -pe novaall 32-
export CALCDIR="/usr/calc/${JOB_ID}/"
export JOBDIR="${PWD}/"
export OMP_NUM_THREADS=$NUM_SLOTS
mkdir $CALCDIR
rsync -av $JOBDIR  $CALCDIR
cd $CALCDIR
phono3py --fc3 --fc2 --dim="1 1 1" --dim_fc2="2 2 2" --mesh="12 12 12" --nac -c POSCAR --br --write_gamma --gp="GP" --ts=300 --pa="0 1/2 1/2 1/2 0 1/2 1/2 1/2 0"
cd $JOBDIR
rsync -av $CALCDIR  $JOBDIR
rm -rf $CALCDIR

