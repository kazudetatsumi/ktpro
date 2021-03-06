#!/usr/bin/env python
# This object reads an INS 4D data and mask h5py file
# and show 2D image of theese projections.
# This script is written to create a neat figure of the additional mask
# investigation.
# Kazuyoshi TATSUMI 2021/03/06
import numpy as np
import sys
import matplotlib.pyplot as plt
import h5py


class Projection:

    def __init__(self, datafile, maskfile, wholerange):
        self.datafile = datafile
        self.maskfile = maskfile
        self.wholerange = wholerange

    def process(self):
        f = h5py.File(self.datafile, 'r')
        g = h5py.File(self.maskfile, 'r')
        self.data4_org = np.array(f["data4"])
        self.condition = np.array(g["condition"])
        self.data4 = self.data4_org*self.condition
        self.plot_projectedintensity()

    def plotter(self, lx, ly, data, vn, hn, cn):
        ax = self.fig.add_subplot(vn,  hn, cn)
        ax.pcolor(np.transpose(data), vmax=np.max(data), cmap='gray')
        ax.set_xlabel(lx)
        ax.set_ylabel(ly)
        ax.xaxis.set_label_coords(0.5, 1.145)
        ax.tick_params(direction="in", color="white", top=True, labeltop=True,
                       labelbottom=False)
        ax.axis('tight')

    def plot_projectedintensity(self):
        self.fig = plt.figure(figsize=(12, 9))
        self.fig.suptitle("projected 4D INS data and unmasked elements from "
                          + self.maskfile)
        self.plotter('qx', 'qy', np.sum(self.data4_org[:, :, :, 0], axis=2),
                     3, 3, 1)
        self.plotter('qx', 'qy', np.sum(self.data4[:, :, :, 0], axis=2),
                     3, 3, 2)
        self.plotter('qx', 'qy', np.sum(np.sum(self.condition, axis=2), axis=2),
                     3, 3, 3)
        self.plotter('qy', 'E', np.sum(self.data4_org[:, :, 0, :], axis=0),
                     3, 3, 4)
        self.plotter('qy', 'E', np.sum(self.data4[:, :, 0, :], axis=0),
                     3, 3, 5)
        self.plotter('qy', 'E', np.sum(np.sum(self.condition, axis=0), axis=1),
                     3, 3, 6)
        self.plotter('qx', 'E', np.sum(self.data4_org[:, :, 0, :], axis=1),
                     3, 3, 7)
        self.plotter('qx', 'E', np.sum(self.data4[:, :, 0, :], axis=1),
                     3, 3, 8)
        self.plotter('qx', 'E', np.sum(np.sum(self.condition, axis=1), axis=1),
                     3, 3, 9)


def samplerun():
    datafile = "/home/kazu/desktop/200522/Ei24/fineq/26m/eliminated_data.hdf5"
    maskfile = "/home/kazu/desktop/200522/Ei24/fineq/add_random_mask/26m/" +\
               "try595_26m/condition.hdf5"
    prj = Projection(datafile, maskfile, 'dummy')
    prj.process()
    plt.show()


samplerun()
