#!/usr/bin/env python
# This script tries to form a densit esimation of qens
# using the bin width distribution of the adaptive KDE
# and the usual qens histogram wiht several corections.
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/kazu/ktpro")
from qens_kde_results_odata_divided_by_idata_class\
    import odata_divided_by_idata as odbi
from get_qlist_nova_class import get_qlist as gq
pwd = os.getcwd()


def getsspectra(sfile, qmin, qmax):
    prj = gq(pklfile=sfile)
    prj.read_pkl()
    prj.spect(qmin, qmax)
    return prj.ene, prj.spectra


def getdivspectra(divfile):
    prj = gq(pklfile=divfile)
    prj.read_pkl()
    return(prj.dataset['xo'], prj.dataset['yr_ssvk'])


def getbandwidth(kf):
    proj = odbi(kf, 'dummy')
    return proj.read_pkl(kf)['y_ssvk']


def balloon(ky, sy):
    bw = np.interp(sy[0], ky[1], ky[2])
    idx = sy[1].nonzero()
    sy_nz = sy[1][idx]
    t_nz = sy[0][idx]
    yv = np.zeros_like(sy[0])
    dt = min(np.diff(sy[0]))
    for k in range(yv.shape[0]):
        yv[k] = np.sum(sy_nz * dt * Gauss(sy[0][k]-t_nz, bw[k]))
    yv = yv * np.sum(ky[0]) / np.sum(yv)
    return yv


def Gauss(x, w):
    y = 1 / (2 * np.pi)**2 / w * np.exp(-x**2 / 2 / w**2)
    return y


def testrun():
    qrange = pwd.split('q')[1]
    qmin = float(qrange.split('-')[0])
    qmax = float(qrange.split('-')[1])
    sprefix = "/home/kazu/desktop/210108/Tatsumi/srlz/0000025/"
    allqf = sprefix + "run6202s.pkl"
    kf = "./qens_run6202united_kde_results_on_data_qsel.pkl"
    kdivf = "./qens_kde_o_divided_by_i_6202.pkl"
    elim = [-0.03, 0.07]
    ky = getbandwidth(kf)
    sy = getsspectra(allqf, qmin, qmax)
    syb = balloon(ky, sy)
    divy = getdivspectra(kdivf)
    smask = np.where((sy[0] > elim[0]) & (sy[0] <= elim[1]))
    kmask = np.where((ky[1] > elim[0]) & (ky[1] <= elim[1]))
    dmask = np.where((divy[0] > elim[0]) & (divy[0] <= elim[1]))
    plt.plot(sy[0][smask], syb[smask])
    plt.plot(ky[1][kmask], ky[0][kmask])
    plt.plot(divy[0][dmask], divy[1][dmask]*10.0)
    plt.yscale('log')
    plt.show()



testrun()
