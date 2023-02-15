#!/usr/bin/env python
# This script tries to form a density esimation of qens
# using the bin width distribution of the adaptive KDE
# and the usual qens histogram wiht several corections.
import os
import numpy as np
import pickle
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/kazu/ktpro")
from qens_kde_results_odata_divided_by_idata_class\
    import odata_divided_by_idata as odbi
from get_qlist_nova_class import get_qlist as gq
pwd = os.getcwd()
from qens_fit_class import qens_fit as qf
###elim = [-0.03, 0.07]
elim = [-0.08, 0.20]
###import matplotlib
###matplotlib.use('Agg')


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
    yv = yv * np.max(ky[0]) / np.max(yv)
    return yv


def Gauss(x, w):
    y = 1 / (2 * np.pi)**2 / w * np.exp(-x**2 / 2 / w**2)
    return y


def plotter(sy, syb, ky, divy, elim, runno, fig):
    smask = np.where((sy[0] > elim[0]) & (sy[0] <= elim[1]))
    kmask = np.where((ky[1] > elim[0]) & (ky[1] <= elim[1]))
    dmask = np.where((divy[0] > elim[0]) & (divy[0] <= elim[1]))
    plt.plot(sy[0][smask], sy[1][smask], label='chk sdata')
    ####plt.plot(sy[0][smask], syb[smask], label='balloon estimate')
    ####plt.plot(ky[1][kmask], ky[0][kmask], label='kde on outgoing neutrons')
    ####plt.plot(divy[0][dmask], divy[1][dmask]/np.max(divy[1][dmask])*np.max(syb[smask]), label='kde divided')
    ####plt.yscale('log')
    ####plt.legend()
    #plt.savefig('balloon_run'+runno+'.png')
    #plt.show()
    plt.show()


def save_pkl(sy, syb, pklfile):
    dataset = {}
    dataset['sy'] = sy
    dataset['syb'] = syb
    with open(pklfile, 'wb') as f:
        pickle.dump(dataset, f, -1)


def eachrunno_read_pkl(runno):
    pklfile = "./qens_run" + runno + "_balloon.pkl"
    with open(pklfile, 'rb') as f:
        dataset = pickle.load(f)
    return dataset['sy'][0], dataset['syb']


def eachrunno(runno, fig):
    ###qrange = pwd.split('q')[1]
    ###print(qrange)
    ###qmin = float(qrange.split('-')[0])
    ###qmax = float(qrange.split('-')[1])
    qmin = 0.55
    qmax = 0.70
    ###sprefix = "/home/kazu/desktop/210108/Tatsumi/srlz/0000001/"
    sprefix = "/home/kazu/desktop/210108/Tatsumi/srlz/0000025/"
    allqf = sprefix + "run" + runno + "s.pkl"
    kf = "./qens_run" + runno + "united_kde_results_on_data_qsel.pkl"
    kdivf = "./qens_kde_o_divided_by_i_" + runno + ".pkl"
    pklfile = "./qens_run" + runno + "_balloon.pkl"
    print(allqf)
    print(kf)
    ky = getbandwidth(kf)
    sy = getsspectra(allqf, qmin, qmax)
    syb = balloon(ky, sy)
    divy = getdivspectra(kdivf)
    plotter(sy, syb, ky, divy, elim, runno, fig)
    ###save_pkl(sy, syb, pklfile)
    return sy[0], syb


def run():
    ###if len(sys.argv) >= 2:
    ###    runnot = sys.argv[1]
    runnot = "6202"
    fig = plt.figure()
    xt, yt = eachrunno(runnot, fig)
    xd, yd = eachrunno("6204", fig)
    ###plt.savefig('balloon_run.png')
    proj = qf('dummy', 'dummy', elim, showplot=False, leastsq=False)
    proj.icorr()
    proj.x_tf, proj.y_tf = proj.limit2(xt, yt, elim)
    proj.x_df, proj.y_df = proj.limit2(xd, yd, elim)
    proj.correction()
    proj.bg = 0.
    fig = plt.figure()
    #proj.optimize(variables=[0.8, 0.01, 0.24, 0.0002, 0.001, 1.2], figname='balloon_fit.png')
    ###proj.optimize(variables=[0.8, 0.01, 0.24, 0.0002], figname='balloon_fit_4params.png')


def fitrun():
    if len(sys.argv) >= 2:
        runnot = sys.argv[1]
    xt, yt = eachrunno_read_pkl(runnot)
    xd, yd = eachrunno_read_pkl("6204")
    proj = qf('dummy', 'dummy', elim, showplot=False, leastsq=False)
    proj.icorr()
    proj.x_tf, proj.y_tf = proj.limit2(xt, yt, elim)
    proj.x_df, proj.y_df = proj.limit2(xd, yd, elim)
    proj.correction()
    proj.bg = proj.y_tf[-1]
    fig = plt.figure()
    #proj.optimize(variables=[0.387, 0.0250, 0.132, 0.00407, 0.437, 0.628], figname='balloon_fit.png')
    #proj.optimize(variables=[0.655, 0.0129, 0.200, 0.00208, 0.171], figname='balloon_fit.png')
    proj.optimize(variables=[0.655, 0.0129, 0.200, 0.00208], figname='balloon_fit.png')
    #proj.optimize(variables=[0.8, 0.0005, 0.2, 0.0002], figname='balloon_fit.png')


run()
