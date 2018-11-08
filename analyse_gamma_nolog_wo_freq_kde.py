#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py
import math
from scipy import stats

homedir = "/home/kazu/"
cdir = homedir + "/asi3n4/phono3py_112_fc2_334_sym_monk_shift/"
sdir = homedir + "/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/"
gdir = homedir + "/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/"
Temp = 300
max_freq = 5
fs = 9
#c = cdir + "noiso/kappa-m141416.noiso.hdf5"
c = cdir + "noiso/kappa-m101014.noiso.hdf5"
#s = sdir + "noiso/kappa-m141432.noiso.hdf5"
s = sdir + "noiso/kappa-m101026.noiso.hdf5"
g = gdir + "noiso/kappa-m181818.hdf5"
#cv = cdir + "qpoints_m141416.hdf5"
cv = cdir + "qpoints.hdf5"
sv = sdir + "qpoints.hdf5"
gv = gdir + "qpoints.hdf5"


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
    # print np.array(omega_k).shape
    # print np.array(gamma_k).shape
    # print np.array(qpoint_k).shape
    # print omega_k[0,0:32]

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
                if k % 3 == 0 or k % 3 == 1:
                    amp += abs(edata[i, k, j])**2
            ampdata[i, j] = amp
    return(ampdata)


def trans(qpoint, edata):
    datasize = edata.shape
    longidata = np.zeros((datasize[0], datasize[2]))
    for i in range(0, datasize[0]):
        norm = (qpoint[i, 0]**2 + qpoint[i, 1]**2 + qpoint[i, 2]**2)**0.5
        if norm > 0:
            for j in range(0, datasize[2]):
                longi = 0
                for k in range(0, datasize[1]):
                    if k % 3 == 0:
                        longi += (abs(edata[i, k, j]*qpoint[i, 0] +
                                  edata[i, k + 1, j] * qpoint[i, 1] +
                                  edata[i, k + 2, j] * qpoint[i, 2])/norm)**2
                longidata[i, j] = longi
    return(longidata)


def kde(m1, m2, n, xlabel, ylabel):
    print m1.shape
    print m2.shape
    xmin = m1.min()
    xmax = m1.max()
    ymin = m2.min()
    ymax = m2.max()
    ymin = 0
    ymax = 200
    #xmin = 0
    #xmax = 35
    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([m1, m2])
    kernel = stats.gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)
    print Z.max()
    print Z.min()
    plt.subplot(2, 2, n)
    plt.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r, extent=[xmin, xmax, ymin, ymax], aspect='auto')
    plt.plot(m1, m2, 'k.', markersize=1)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)


def select_mode(omega, gamma, sqamp, lamp):
    freqs = []
    gammas = []
    sqamps = []
    lamps = []
    for f, g, s, l in zip(omega, gamma, sqamp, lamp):
        condition = f < max_freq
        _f = np.extract(condition, f)
        _g = np.extract(condition, g)
        _s = np.extract(condition, s)
        _l = np.extract(condition, l)
        freqs += list(_f)
        gammas += list(_g)
        sqamps += list(_s)
        lamps += list(_l)
    f1 = np.array(freqs).ravel()
    g1 = np.array(gammas).ravel()
    s1 = np.array(sqamps).ravel()
    l1 = np.array(lamps).ravel()
    return(f1, g1, s1, l1)


def caserun(casefile, vcasefile, n, phase):
    omegak, gammak, qpointk = parse_gamma(casefile, Temp)
    omegaq, eigenvecq, qpointq = parse_eigenvec(vcasefile)
    check(omegak, omegaq)
    check(qpointk, qpointq)
    sqamp = project_eigenvec(eigenvecq)
    longiamp = trans(qpointq, eigenvecq)
    omega1d, gamma1d, sqamp1d, longiamp1d = select_mode(omegak, gammak, sqamp,
                                                        longiamp)
    kde(longiamp1d[3:], 1.0 / (4 * np.pi * gamma1d[3:]), n,'fraction of longitudinal components', 'tau')
    kde(sqamp1d[3:], 1.0 / (4 * np.pi * gamma1d[3:]), n+2, 'fraction of eigenvector components on x-y plane', 'tau')
    #kde(omega1d[3:], 1.0 / (4 * np.pi * gamma1d[3:]), n, 'omega', 'tau')
    #plt.subplot(2, 2, n)
    #plt.scatter(longiamp1d, 1.0 / (4 * np.pi * gamma1d),
    #            linewidth=0.01, s=4, color='k', marker='o',
    #            label=phase)
    #plt.tick_params(which='both', tickdir='in')
    #plt.ylim(20, 110)
    #plt.yticks([20, 20, 40, 60, 80, 100])
    #plt.xlim(-0.01, 1.01)
    #plt.xlabel('fraction of longitudinal components')
    #plt.ylabel('tau')
    #plt.subplot(2, 2, n+2)
    #plt.scatter(sqamp1d, 1.0 / (4 * np.pi * gamma1d),
    #            linewidth=0.01, s=4, color='k', marker='o', label=phase)
    #plt.tick_params(which='both', tickdir='in')
    #n_left = 0
    #n_right = 0
    #for s in sqamp1d:
    #    if s > 0.5:
    #        n_right += 1
    #    else:
    #        n_left += 1
    #print "num of left states:", n_left
    #print "num of right states:", n_right

    #plt.ylim(20, 110)
    #plt.yticks([20, 40, 60, 80, 100])
    #plt.xlim(-0.01, 1.01)
    #plt.xlabel('fraction of eigenvector components on x-y plane')
    #plt.ylabel('tau')


def run():
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.figure(figsize=(9, 9))
    caserun(c, cv, 1, "alpha")
    caserun(s, sv, 2, "beta")
    #caserun(g, gv, 3, "gamma")
    #plt.savefig("tst_plot.eps")

run()
plt.show()
