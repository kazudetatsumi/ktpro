#!/usr/bin/env python
import numpy as np
import h5py
import random
import os
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
psize = comm.Get_size()


class Add_Random_Mask:

    def __init__(self, timelist, num_try, workdir, tail):
        self.timelist = timelist
        self.num_try = num_try
        self.workdir = workdir
        self.tail = tail

    def process(self):
        for time in self.timelist:
            self.outfile = self.workdir + "/" + str(time) + self.tail
            f = h5py. File(self.outfile, 'r')
            data4_org = np.array(f["data4"])
            condition_org = np.array(f["condition"])
            ylen = data4_org.shape[1]
            elen = data4_org.shape[3]
            mrank = psize - rank - 1
            if self.num_try >= mrank + 1:
                for tryidx in range(mrank + 0,
                                    self.num_try -
                                    ((self.num_try-mrank-1) % psize) +
                                    psize - 1,
                                    psize):
                    condition = np.array(f["condition"])
                    dirname = "./try" + str(tryidx) + "_" + str(time) +\
                              self.tail.split("/")[0]
                    savefile = "./try" + str(tryidx) + "_" + str(time) +\
                               self.tail
                    os.system("mkdir " + dirname)
                    yc = random.random()*ylen
                    ec = random.random()*elen
                    radius = (ylen**2 + elen**2)**0.5*0.5*random.random()
                    for yindx in range(0, condition.shape[1]):
                        for eindx in range(0, condition.shape[3]):
                            if abs(yindx - yc)**2 + abs(eindx - ec)**2 \
                               < radius**2:
                                condition[:, yindx, :, eindx] = False
                    data4 = data4_org*condition
                    num_ele = np.sum(condition)
                    num_ele_org = np.sum(condition_org)
                    frac_add_mask = (num_ele_org - num_ele)*1.0 / num_ele_org
                    self.save_eliminated_data_hdf5(savefile, data4, condition,
                                                   frac_add_mask)

    def save_eliminated_data_hdf5(self, savefile, data, condition,
                                  frac_add_mask):
        with h5py.File(savefile, 'w') as hf:
            hf.create_dataset('data4', data=data)
            hf.create_dataset('condition', data=condition)
            hf.create_dataset('frac_add_mask', data=frac_add_mask)


def sample_run():
    timelist = [26]
    num_try = 10
    workdir = "/home/kazu/desktop/200522/Ei24/fineq/condparam07/"
    tail = "m/eliminated_data.hdf5"
    projectset = Add_Random_Mask(timelist, num_try, workdir, tail)
    projectset.process()


sample_run()
