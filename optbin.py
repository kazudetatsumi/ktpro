#!/usr/bin/env python
# sample program of using optbinwidth4d_wholefort_mpi_class
# u should modifie datafile and condfile 
# mpirun -np num_proc optbin.py
import sys
sys.path.append("/home/kazu/ktpro")
import optbinwidth4d_wholefort_mpi_class as ob


def run():
    datafile = "/home/kazu/desktop/200204/fine/hourbyhour/" +\
                "ortho_opt_without_mask/10h/eliminated_data.hdf5"
    condfile = "./condition.hdf5"
    usecond = False
    projectset = ob.Opt_Bin_Width(datafile, condfile, usecond)
    data, condition, maxw = projectset.preprocess()
    data = data*condition*1.0
    projectset.calc_cost4d_f90(data, condition, maxw)


run()
