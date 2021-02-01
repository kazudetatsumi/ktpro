#!/usr/bin/env python
# This is an exercise for writing a class object.
# This object reads an INS 4D data and mask h5py file
# and adds a cylindrical mask whose position and diameter are set at random.
# For preparation of input files for investigation of  the effect of
# artificial mask on the optimized bin-widths.
# This object is adaptive to mpi parallelism.
# changed not to save INS intensities but to save mask info more.
# Kazuyoshi TATSUMI 2021/01/26


import numpy as np
import h5py
import random
import os
from mpi4py import MPI
from matplotlib import pyplot as plt


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
psize = comm.Get_size()


class Add_Random_Mask:

    def __init__(self, timelist, num_try, workdir, tail, maxradius):
        self.timelist = timelist
        self.num_try = num_try
        self.workdir = workdir
        self.tail = tail
        self.maxradius = maxradius

    def process(self, axis=1, itry=0):
        for time in self.timelist:
            outfile = self.workdir + "/" + str(time) + self.tail
            f = h5py. File(outfile, 'r')
            data4_org = np.array(f["data4"])
            condition_org = np.array(f["condition"])
            ylen = data4_org.shape[axis]
            elen = data4_org.shape[3]
            mrank = psize - rank - 1
            if self.num_try >= mrank + 1:
                for tryidx in range(itry + mrank + 0,
                                    itry + self.num_try -
                                    ((self.num_try-mrank-1) % psize) +
                                    psize - 1,
                                    psize):
                    condition = np.array(f["condition"])
                    dirname = "./try" + str(tryidx) + "_" + str(time) +\
                              self.tail.split("/")[0]
                    savefile = dirname + "/condition.hdf5"
                    figfile = "data4_" + str(time) + "_" + str(tryidx) + ".png"
                    os.system("mkdir " + dirname)
                    radius = self.maxradius*random.random()
                    yc = radius + random.random()*(ylen - 2*radius)
                    ec = radius + random.random()*(elen - 2*radius)
                    for yindx in range(0, condition.shape[axis]):
                        for eindx in range(0, condition.shape[3]):
                            if abs(yindx - yc)**2 + abs(eindx - ec)**2 \
                               < radius**2:
                                if axis == 0:
                                    condition[yindx, :, :, eindx] = False
                                elif axis == 1:
                                    condition[:, yindx, :, eindx] = False
                                elif axis == 2:
                                    condition[:, :, yindx, eindx] = False
                    data4 = data4_org*condition
                    num_ele = np.sum(condition)
                    num_ele_org = np.sum(condition_org)
                    frac_add_mask = (num_ele_org - num_ele)*1.0 / num_ele_org
                    self.save_condition_hdf5(savefile, data4, condition,
                                             frac_add_mask, yc, ec, radius)
                    self.plot_projectedintensity(condition, data4)
                    plt.savefig(figfile)

    def save_condition_hdf5(self, savefile, data, condition, frac_add_mask,
                            yc, ec, radius):
        with h5py.File(savefile, 'w') as hf:
            hf.create_dataset('condition', data=condition)
            hf.create_dataset('frac_add_mask', data=frac_add_mask)
            hf.create_dataset('yc', data=yc)
            hf.create_dataset('ec', data=ec)
            hf.create_dataset('radius', data=radius)

    def plotter(self, lx, ly, data, vn, hn, cn):
        ax = self.fig.add_subplot(vn,  hn, cn)
        ax.pcolor(np.transpose(data), vmax=np.max(data), cmap='jet')
        ax.set_xlabel(lx)
        ax.set_ylabel(ly)
        ax.xaxis.set_label_coords(0.5, 1.145)
        ax.tick_params(direction="in", color="white", top=True, labeltop=True,
                       labelbottom=False)
        ax.axis('tight')

    def plot_projectedintensity(self, condition, data4):
        self.fig = plt.figure(figsize=(9, 9))
        self.fig.suptitle("crosssections of 4D INS data from" + self.workdir)
        self.plotter('qx', 'qy', np.sum(data4[:, :, :, 0], axis=2), 3, 2, 1)
        self.plotter('qx', 'qy', np.sum(np.sum(condition, axis=2), axis=2),
                     3, 2, 2)
        self.plotter('qy', 'E', np.sum(data4[:, :, 0, :], axis=0), 3, 2, 3)
        self.plotter('qy', 'E', np.sum(np.sum(condition, axis=0), axis=1),
                     3, 2, 4)
        self.plotter('qx', 'E', np.sum(data4[:, :, 0, :], axis=1), 3, 2, 5)
        self.plotter('qx', 'E', np.sum(np.sum(condition, axis=1), axis=1),
                     3, 2, 6)


#def sample_run():
#    timelist = [26]
#    num_try = 10
#    workdir = "/home/kazu/desktop/200522/Ei24/fineq/condparam07/"
#    tail = "m/eliminated_data.hdf5"
#    maxradius = 10.0 # per step
#    projectset = Add_Random_Mask(timelist, num_try, workdir, tail, maxradius)
#    projectset.process(axis=1)
#
#
#sample_run()
