#!/usr/bin/env python
# This object reads an INS 4D data and mask h5py file
# and show 2D image of theese projections.
# This script is written to create a neat figure of the additional mask
# investigation.
# Kazuyoshi TATSUMI 2021/03/06
import numpy as np
import sys
sys.path.append('/home/kazu/ktpro/')
import crosssection_of_4ddata_class as c4c
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py
from matplotlib.colors import LinearSegmentedColormap


class SProjection(c4c.CROSS):

    def __init__(self, datafile, maskfile, cross=False, orthotope_lims=None,
                 wholeranges=None, devs=np.array([1, 1, 1,]), cpos=None, hvlofs=False,
                 binwidths=None, common_lims=None, filled=False):
        self.datafile = datafile
        self.maskfile = maskfile
        self.cross = cross
        self.orthotope_lims = orthotope_lims
        self.wholeranges = wholeranges
        self.devs = devs
        self.cpos = cpos
        self.hvlofs = hvlofs
        self.binwidths = binwidths
        self.common_lims = common_lims
        self.filled = filled

    def process(self):
        f = h5py.File(self.datafile, 'r')
        g = h5py.File(self.maskfile, 'r')
        self.data4_org = np.array(f["data4"])
        self.condition_org = np.array(f["condition"])
        self.condition = np.array(g["condition"])
        self.data4 = self.data4_org*self.condition
        if self.cross:
            self.data4[self.condition <= 0.1] = -10000.0
            super(SProjection, self).create_fig()
            super(SProjection, self).plot_crosssection(3)
            self.data4_org[self.condition_org <= 0.1] = -10000.0
            self.data4 = self.data4_org
            super(SProjection, self).plot_crosssection(2)
        else:
            super(SProjection, self).create_fig()
            self.plot_projectedintensity()

    def create_cmap(self):
        self.custom_cmap = cm.get_cmap('jet')
        self.custom_cmap.set_under(color='black')

    def prjplotter(self, lx, ly, data, vn, hn, cn):
        self.create_cmap()
        ax = self.fig.add_subplot(vn,  hn, cn)
        ax.pcolor(np.transpose(data), vmin=0.0, vmax=np.max(data), cmap=self.custom_cmap)
        ax.set_xlabel(lx)
        ax.set_ylabel(ly)
        ax.xaxis.set_label_coords(0.5, 1.145)
        ax.tick_params(direction="in", color="white", top=True, labeltop=True,
                       labelbottom=False)
        ax.axis('tight')


    def add_negative(self, axis1, axis2):
        self.sumcondition = np.sum(np.sum(self.condition, axis=axis1), axis=axis2)
        self.sumcondition_org = np.sum(np.sum(self.condition_org, axis=axis1), axis=axis2)
        self.sumdata4_org = np.sum(np.sum(self.data4_org, axis=axis1), axis=axis2)
        self.sumdata4 = np.sum(np.sum(self.data4, axis=axis1), axis=axis2)
        self.sumdata4_org[self.sumcondition_org == 0] = -100.0
        self.sumdata4[self.sumcondition == 0] = -100.0
        self.sumcondition = self.sumcondition - 1

    def plot_projectedintensity(self):
        #self.fig = plt.figure(figsize=(12, 9))
        #self.fig.suptitle("projected 4D INS data and unmasked elements from "
        #                  + self.maskfile)
        lims_int = super(SProjection, self).lims_int(self.orthotope_lims)
        cposinfo = ""

        self.add_negative(2, 2)
        alpha = np.max(self.sumdata4)/np.max(self.sumdata4_org)
        self.plotter(lims_int[0, :], lims_int[1, :], '$q_a(rlu)$',
                     '$q_b(rlu)$', self.sumdata4_org,
                     4, 3, 1, self.devs[0], self.wholeranges[0, :], 
                     self.wholeranges[1, :], cposinfo)
        self.plotter(lims_int[0, :], lims_int[1, :], '$q_a(rlu)$',
                     '$q_b(rlu)$', self.sumdata4,
                     4, 3, 2, alpha, self.wholeranges[0, :], 
                     self.wholeranges[1, :], cposinfo)
        self.plotter(lims_int[0, :], lims_int[1, :], '$q_a(rlu)$',
                     '$q_b(rlu)$', self.sumcondition,
                     4, 3, 3, self.devs[0], self.wholeranges[0, :], 
                     self.wholeranges[1, :], cposinfo)

        self.add_negative(0, 1)
        alpha = np.max(self.sumdata4)/np.max(self.sumdata4_org)
        self.plotter(lims_int[1, :], lims_int[3, :], '$q_b(rlu)$',
                     'E(meV)', self.sumdata4_org,
                     4, 3, 4, self.devs[1], self.wholeranges[1, :], 
                     self.wholeranges[3, :], cposinfo)
        self.plotter(lims_int[1, :], lims_int[3, :], '$q_b(rlu)$',
                     'E(meV)', self.sumdata4,
                     4, 3, 5, alpha, self.wholeranges[1, :], 
                     self.wholeranges[3, :], cposinfo)
        self.plotter(lims_int[1, :], lims_int[3, :], '$q_b(rlu)$',
                     'E(meV)', self.sumcondition,
                     4, 3, 6, self.devs[1], self.wholeranges[1, :], 
                     self.wholeranges[3, :], cposinfo)

        self.add_negative(1, 1)
        alpha = np.max(self.sumdata4)/np.max(self.sumdata4_org)
        self.plotter(lims_int[0, :], lims_int[3, :], '$q_a(rlu)$',
                     'E(meV)', self.sumdata4_org,
                     4, 3, 7, self.devs[2], self.wholeranges[0, :], 
                     self.wholeranges[3, :], cposinfo)
        self.plotter(lims_int[0, :], lims_int[3, :], '$q_a(rlu)$',
                     'E(meV)', self.sumdata4,
                     4, 3, 8, alpha, self.wholeranges[0, :], 
                     self.wholeranges[3, :], cposinfo)
        self.plotter(lims_int[0, :], lims_int[3, :], '$q_a(rlu)$',
                     'E(meV)', self.sumcondition,
                     4, 3, 9, self.devs[2], self.wholeranges[0, :], 
                     self.wholeranges[3, :], cposinfo)

        self.add_negative(0, 0)
        self.plotter(lims_int[2, :], lims_int[3, :], '$q_c(rlu)$',
                     'E(meV)', self.sumdata4_org,
                     4, 3, 10, self.devs[2], self.wholeranges[2, :], 
                     self.wholeranges[3, :], cposinfo)
        self.plotter(lims_int[2, :], lims_int[3, :], '$q_c(rlu)$',
                     'E(meV)', self.sumdata4,
                     4, 3, 11, self.devs[2], self.wholeranges[2, :], 
                     self.wholeranges[3, :], cposinfo)
        self.plotter(lims_int[2, :], lims_int[3, :], '$q_c(rlu)$',
                     'E(meV)', self.sumcondition,
                     4, 3, 12, self.devs[2], self.wholeranges[2, :], 
                     self.wholeranges[3, :], cposinfo)


def samplerun():
    datafile = "/home/kazu/desktop/200522/Ei24/fineq/26m/eliminated_data.hdf5"
    maskfile = "/home/kazu/desktop/200522/Ei24/fineq/add_random_mask/26m/" +\
               "try138_26m/condition.hdf5"
    prj = SProjection(datafile, maskfile, cross=False)
    prj.orthotope_lims = np.array([[0, 228], [0, 202], [0, 9], [0, 139]])*1.0
    orgbinwidths = np.array([2, 3, 2, 4])
    prj.binwidths = np.array([1, 1, 1, 1])*1.
    prj.wholeranges = np.array([[0.01, 2.31], [-0.67, 1.35], [-0.16, 0.18],
                                [10.0, 21.20]])
    prj.mb = 0.1
    prj.showorthob = False
    prj.cpos = np.array([75, 39, 2, 6])*orgbinwidths

    prj.process()
    plt.savefig("example_mask.png")
    plt.show()


samplerun()
