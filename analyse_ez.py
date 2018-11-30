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
fig, ax = plt.subplots(2, 6, figsize=(24, 8))


def parse_gamma(filename, temp):
    f = h5py.File(filename, 'r')
    temperature = f["temperature"].value
    i = 0
    for t in temperature:
        if t == temp:
            tindex = i
        i += 1
    gamma_k = f["gamma"][tindex, ]
    omega_k = f["frequency"][:, :]
    qpoint_k = f["qpoint"][:]
    return(omega_k, gamma_k, qpoint_k)


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


def parse_kpoints(qpfile, M):
    f = h5py.File(qpfile)
    qpt = f["qpoint"]
    kpt = np.dot(qpt, np.transpose(np.linalg.inv((M))))
    return(kpt)


def band_by_band(omega, gamma, sqamp, kcoordx, kcoordz, phase, n, m,
                 markstyle):
    freqs = []
    gammas = []
    sqamps = []
    kcdxs = []
    kcdzs = []
    for f, g, s, kx, kz in zip(omega, gamma, sqamp, kcoordx, kcoordz):
        condition = f < max_freq
        _f = np.extract(condition, f)
        _g = np.extract(condition, g)
        _s = np.extract(condition, s)
        _kx = np.extract(condition, kx)
        _kz = np.extract(condition, kz)
        freqs += list(_f)
        gammas += list(_g)
        sqamps += list(_s)
        kcdxs += list(_kx)
        kcdzs += list(_kz)
    f1 = np.array(freqs).ravel()
    g1 = np.array(gammas).ravel()
    s1 = np.array(sqamps).ravel()
    kx1 = np.array(kcdxs).ravel()
    kz1 = np.array(kcdzs).ravel()
    sc = ax[n, m*2].scatter(f1, kx1, c=s1, linewidth=0.01,
                            s=18, label=phase, cmap='gnuplot',
                            marker=markstyle)
    ax[n, m*2].set_ylim(-0.10, 0.10)
    ax[n, m*2].set_title(phase)
    ax[n, m*2].set_xlabel('omega / THz')
    ax[n, m*2].set_ylabel('kx / Angs-1')
    sc = ax[n, 1 + m*2].scatter(f1, kz1, c=s1, linewidth=0.01,
                                s=18, label=phase, cmap='gnuplot',
                                marker=markstyle)
    ax[n, 1 + m*2].set_ylim(-0.10, 0.10)
    ax[n, 1 + m*2].set_title(phase)
    ax[n, 1 + m*2].set_ylabel('kz / Angs-1')
    if markstyle == "o":
        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        fig.colorbar(sc, cax=cbar_ax)


def caserun(casefile, vcasefile, A, n, phase):
    omegak, gammak, qpointk = parse_gamma(casefile, Temp)
    omegaq, eigenvecq, qpointq = parse_eigenvec(vcasefile)
    check(omegak, omegaq)
    check(qpointk, qpointq)
    kpt = parse_kpoints(vcasefile, A)
    sqamp = project_eigenvec(eigenvecq)
    band_by_band(omegak[:, 0], gammak[:, 0], sqamp[:, 0], kpt[:, 0],
                 kpt[:, 2], phase, n, 0, "o")
    band_by_band(omegak[:, 1], gammak[:, 1], sqamp[:, 1], kpt[:, 0],
                 kpt[:, 2], phase, n, 1, "v")
    band_by_band(omegak[:, 2], gammak[:, 2], sqamp[:, 2], kpt[:, 0],
                 kpt[:, 2], phase, n, 2, "d")


def run():
    cA = np.array([[7.8079969443386954, 0.0000000000000000, 0.0000000000000001], [-3.9039984721693477, 6.7619237064685818, 0], [0, 0, 5.6590777347249741]])
    sA = np.array([[7.6595137795552795, 0.0000000000000000, 0.0000000000000001], [-3.8297568897776397, 6.6333335137318326, 0], [0, 0, 2.9247116510287272]])
    caserun(c, cv, cA, 0, "alpha")
    caserun(s, sv, sA, 1, "beta")
    plt.savefig("tst_plot.pdf")

run()



plt.show()
