#!/usr/bin/env python
import sys
sys.path.append("/home/kazu/ktpro")
import numpy as np
import matplotlib.pyplot as plt
from get_resampled_data_mpi_class import Sget_qlist as sgq
from qens_fit_class_kder import runkdenoidata as rkn


class qens_kde_resampled(sgq, rkn):
    def __init__(self, pklfile):
        self.pklfile = pklfile
        self.load_pkl()
        #self.kde(self.spectrab[inb, 0, :], self.spectrab[inb, 2, :])


def testrun():
    inb = 0
    pklfile = 'run6202spectrab.pkl'
    proj = qens_kde_resampled(pklfile)
    proj.kde(proj.spectrab[inb, 0, :], proj.spectrab[inb, 2, :], num=800)
    print(len(proj.y[0]))
    print(proj.spectrab[0, 0, :].shape)
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax2 = ax1.twinx()
    ax1.plot(proj.spectrab[0, 0, :], proj.spectrab[0, 2, :]/np.max(proj.spectrab[0, 2, :]))
    ax1.plot(proj.y[1], proj.y[0]/np.max(proj.y[0]))
    ax2.plot(proj.y[1], proj.y[2], c='k')
    plt.savefig('tmp.png')

    

#testrun()
