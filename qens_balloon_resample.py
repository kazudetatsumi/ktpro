#!/usr/bin/env python
# This script balloon-estimates the re-sampled qens 1D intensity profiles with
# the already obtained distributions of kernel band widths, and fits the
# estimated target density with the estimated device density.
# Kazuyoshi TATSUMI 2023/02/17
import os
import numpy as np
import pickle
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/kazu/ktpro")
from qens_kde_results_odata_divided_by_idata_class\
    import odata_divided_by_idata as odbi
#from get_qlist_nova_class import get_qlist as gq
from get_resampled_data_class import Sget_qlist as sgq
pwd = os.getcwd()
from qens_fit_class import qens_fit as qf
elim = [-0.03, 0.07]
import matplotlib
matplotlib.use('Agg')


def getrsspectra(rsfile, inb=0):
    prj = sgq(pklfile=rsfile)
    prj.load_pkl(ispython2=True)
    #return prj.spectrab[0, 0, :], prj.spectrab[0, 1, :]
    return prj.spectrab[inb, 0, :], prj.spectrab[inb, 1, :]


def getdivspectra(divfile):
    prj = sgq(pklfile=divfile)
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
    plt.plot(sy[0][smask], syb[smask], label='balloon estimate')
    plt.plot(ky[1][kmask], ky[0][kmask], label='kde on outgoing neutrons')
    plt.plot(divy[0][dmask], divy[1][dmask]/np.max(divy[1][dmask])*np.max(syb[smask]), label='kde divided')
    plt.yscale('log')
    plt.legend()
    #plt.savefig('balloon_run'+runno+'.png')
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


def eachrunno(runno, fig, inb=0):
    #qrange = pwd.split('q')[1]
    #print(qrange)
    #qmin = float(qrange.split('-')[0])
    #qmax = float(qrange.split('-')[1])
    #sprefix = "/home/kazu/desktop/210108/Tatsumi/srlz/0000001/"
    rsfile = "./run" + runno + "spectrab.pkl"
    kf = "./qens_run" + runno + "united_kde_results_on_data_qsel.pkl"
    #kdivf = "./qens_kde_o_divided_by_i_" + runno + ".pkl"
    #pklfile = "./qens_runr" + runno + "_balloon.pkl"
    #print(rsfile)
    #print(kf)
    ky = getbandwidth(kf)
    sy = getrsspectra(rsfile, inb=inb)
    syb = balloon(ky, sy)
    #divy = getdivspectra(kdivf)
    #plotter(sy, syb, ky, divy, elim, runno, fig)
    #save_pkl(sy, syb, pklfile)
    return sy[0], syb, sy[1]


def run():
    if len(sys.argv) >= 2:
        if sys.argv[1] == 'hist':
            ishist = True
        if sys.argv[1] == 'kde':
            ishist = False
    else:
        ishist = False
    fig = plt.figure()
    Nb = 96
    gammas = np.zeros((Nb, 2))
    np.set_printoptions(threshold=12, linewidth=150, suppress=True)
    print("estimated constants alpha1, gamma1, alpha2, gamma2, delta")
    for inb in range(Nb):
        #print('entering in #', inb)
        xt, yt, yth = eachrunno("6202", fig, inb=inb)
        xd, yd, ydh = eachrunno("6204", fig, inb=inb)
        #plt.savefig('balloon_run.png')
        proj = qf('dummy', 'dummy', elim, showplot=False, leastsq=False, quiet=True)
        proj.icorr()
        if ishist:
            proj.x_tf, proj.y_tf = proj.limit2(xt, yth, elim)
            proj.x_df, proj.y_df = proj.limit2(xd, ydh, elim)
        else:
            proj.x_tf, proj.y_tf = proj.limit2(xt, yt, elim)
            proj.x_df, proj.y_df = proj.limit2(xd, yd, elim)
        proj.correction()
        proj.bg = 0.
        #fig = plt.figure()
        proj.optimize(variables=[0.8, 0.01, 0.24, 0.0002, 0.001, 1.2], figname='balloon_fit_resampled.png')
        gammas[inb, 0] = proj.out[1]
        gammas[inb, 1] = proj.out[3]
    if ishist:
        print('fitting on histogram data')
    else:
        print('fitting on kde data')
    print(np.mean(gammas, axis=0))
    print(np.std(gammas, axis=0))
    #proj.optimize(variables=[0.8, 0.01, 0.24, 0.0002], figname='balloon_fit_4params.png')


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
