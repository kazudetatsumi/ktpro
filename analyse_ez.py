#!/usr/bin/env python
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import h5py
import math
homedir = "/home/kazu/"
cdir = homedir + "/asi3n4/phono3py_112_fc2_334_sym_monk_shift/"
sdir = homedir + "/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/"
gdir = homedir + "/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/"
Temp = 300
max_freq = 5
c = cdir + "noiso/kappa-m101014.noiso.hdf5"
s = sdir + "noiso/kappa-m101026.noiso.hdf5"
g = gdir + "noiso/kappa-m121212.hdf5"
cv = cdir + "qpoints.hdf5"
sv = sdir + "qpoints.hdf5"
gv = gdir + "qpoints.hdf5"

plt.style.use('classic')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['ytick.major.size'] = 5
fig, ax = plt.subplots(2, 8, figsize=(24, 8))


def parse_modekappa(filename, temp):
    f = h5py.File(filename, 'r')
    temperature = f["temperature"].value
    i = 0
    for t in temperature:
        if t == temp:
            tindex = i
        i += 1
    modekappa_k = f["mode_kappa"][tindex, :, :, :]
    omega_k = f["frequency"][:, :]
    qpoint_k = f["qpoint"][:]
    return(omega_k, modekappa_k, qpoint_k)


def parse_eigenvec(filename):
    f = h5py.File(filename, 'r')
    omega_q = f["frequency"][:, :]
    qpoint_q = f["qpoint"][:]
    eigenvec_q = f["eigenvector"][:]
    return(omega_q, eigenvec_q, qpoint_q)


def check(data1, data2):
    data1_1d = np.array(data1).ravel()
    data2_1d = np.array(data2).ravel()
    for d1, d2 in zip(data1_1d, data2_1d):
        if abs(d1 - d2) > 0.1:
            print "Too large difference!"
            print d1, d2
    print "check finished"


def project_eigenvec(edata):
    datasize = edata.shape
    ampdata = np.zeros((datasize[0], datasize[2]))
    for i in range(0, datasize[0]):
        for j in range(0, datasize[2]):
            amp = 0
            for k in range(0, datasize[1]):
                if k % 3 == 2:
                    amp += abs(edata[i, k, j])**2
            ampdata[i, j] = amp
    return(ampdata)


def select_mode(omega, gamma, sqamp):
    freqs = []
    gammas = []
    sqamps = []
    lamps = []
    for f, g, s in zip(omega, gamma, sqamp):
        condition = f < max_freq
        _f = np.extract(condition, f)
        _g = np.extract(condition, g)
        _s = np.extract(condition, s)
        freqs += list(_f)
        gammas += list(_g)
        sqamps += list(_s)
    f1 = np.array(freqs).ravel()
    g1 = np.array(gammas).ravel()
    s1 = np.array(sqamps).ravel()
    return(f1, g1, s1)


def band_by_band(omega, sqamp, mkx, mkz, phase, n, m, markstyle):
    freqs = []
    sqamps = []
    mkxs = []
    mkzs = []
    for f, s, kx, kz in zip(omega, sqamp, mkx, mkz):
        condition = f < max_freq
        _f = np.extract(condition, f)
        _s = np.extract(condition, s)
        _kx = np.extract(condition, kx)
        _kz = np.extract(condition, kz)
        freqs += list(_f)
        sqamps += list(_s)
        mkxs += list(_kx)
        mkzs += list(_kz)
    f1 = np.array(freqs).ravel()
    s1 = np.array(sqamps).ravel()
    mkx1 = np.array(mkxs).ravel()
    mkz1 = np.array(mkzs).ravel()
    sc = ax[n, m*2].scatter(f1, mkx1, c=s1, linewidth=0.01,
                            s=18, label=phase, cmap='gnuplot',
                            marker=markstyle)
    ax[n, m*2].set_ylim(0.0, 1.50)
    ax[n, m*2].set_title(phase)
    ax[n, m*2].set_xlabel('omega / THz')
    ax[n, m*2].set_ylabel('kx / WmK-1')
    sc = ax[n, 1 + m*2].scatter(f1, mkz1, c=s1, linewidth=0.01,
                                s=18, label=phase, cmap='gnuplot',
                                marker=markstyle)
    ax[n, 1 + m*2].set_ylim(0., 1.50)
    ax[n, 1 + m*2].set_title(phase)
    ax[n, 1 + m*2].set_ylabel('kz / WmK-1')
    if markstyle == "o":
        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        fig.colorbar(sc, cax=cbar_ax)


def caserun(casefile, vcasefile, n, phase):
    omegak, modekappak, qpointk = parse_modekappa(casefile, Temp)
    omegaq, eigenvecq, qpointq = parse_eigenvec(vcasefile)
    check(omegak, omegaq)
    check(qpointk, qpointq)
    sqamp = project_eigenvec(eigenvecq)
    band_by_band(omegak[:, 0], sqamp[:, 0], modekappak[:, 0, 0],
                 modekappak[:, 0, 2], phase, n, 0, "o")
    band_by_band(omegak[:, 1], sqamp[:, 1], modekappak[:, 1, 0],
                 modekappak[:, 1, 2], phase, n, 1, "v")
    band_by_band(omegak[:, 2], sqamp[:, 2], modekappak[:, 2, 0],
                 modekappak[:, 2, 2], phase, n, 2, "d")
    rsize = omegak[:, 3:].shape
    romega = np.reshape(omegak[:, 3:], (rsize[0]*rsize[1]))
    rsqamp = np.reshape(sqamp[:, 3:], (rsize[0]*rsize[1]))
    rmodekappakx = np.reshape(modekappak[:, 3:, 0], (rsize[0]*rsize[1]))
    rmodekappakz = np.reshape(modekappak[:, 3:, 2], (rsize[0]*rsize[1]))
    band_by_band(romega, rsqamp, rmodekappakx, rmodekappakz, phase, n, 3, "D")


def run():
    caserun(c, cv, 0, "alpha")
    caserun(s, sv, 1, "beta")
    plt.savefig("tst_plot.pdf")

run()
plt.show()
