#$ -cwd
#$ -V -S /bin/bash
#$ -N phono3py_gp
#$ -pe mpi* 16

phono3py --fc3 --fc2 --dim="1 1 3" --dim_fc2="2 2 4" --mesh="14 14 32" -c POSCAR --nac --br --write_gamma --gp="NUM" --tmax=1600 --tmin=100 --tstep=10 --isotope

