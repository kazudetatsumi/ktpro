#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy import stats
from matplotlib import rc

Temp = 300
nbins = 100
y_max = 0.12
numr = 3
homedir = "/home/kazu/"
cdir = homedir + "asi3n4/phono3py_112_fc2_334_sym_monk_shift/"
sdir = homedir + "bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/"
gdir = homedir + "gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/"
apc = cdir + "gpjob_m101014_fullpp/kappa-m101014.hdf5"
aps = sdir + "gpjob_m101026_fullpp/kappa-m101026.hdf5"
apg = gdir + "gpjob_m181818_fullpp/kappa-m181818.hdf5"
cz = 4*7*3
sz = 2*7*3
gz = 2*7*3


def parse_gamma(filename, temp, max_freq):
    freqs = []
    mode_prop = []
    f = h5py.File(filename, 'r')
    temperature = f["temperature"].value
    i = 0
    for t in temperature:
        if t == temp:
            tindex = i
        i += 1
    gamma = f["gamma"][tindex, ]
    omega = f["frequency"][:, :]

    for freq, g in zip(omega, gamma):
        condition = freq < max_freq
        _freq = np.extract(condition, freq)
        _g = np.extract(condition, g)
        freqs += list(_freq)
        mode_prop += list(_g)
    gamma1 = np.array(mode_prop).ravel()
    omega1 = np.array(freqs).ravel()
    return(omega1, gamma1)


def parse_avepp(filename, max_freq):
    freqs = []
    mode_prop = []
    f = h5py.File(filename, 'r')
    avepp = f["ave_pp"][:, :]
    omega = f["frequency"][:, :]
    for freq, ap in zip(omega, avepp):
        condition = freq < max_freq
        _freq = np.extract(condition, freq)
        _ap = np.extract(condition, ap)
        freqs += list(_freq)
        mode_prop += list(_ap)
    avepp1 = np.array(mode_prop).ravel()
    omega1 = np.array(freqs).ravel()
    return(omega1, avepp1)


def run():
    x = np.array([])
    cc = np.array([])
    cs = np.array([])
    cg = np.array([])
    agc = np.array([])
    ags = np.array([])
    agg = np.array([])
    aac = np.array([])
    aas = np.array([])
    aag = np.array([])
    aagc = np.array([])
    aags = np.array([])
    aagg = np.array([])
    maxf = np.arange(2.0, 35, 0.1)
    for max_freq in maxf:
        omegaapc1, apc1 = parse_avepp(apc, max_freq)
        omegaaps1, aps1 = parse_avepp(aps, max_freq)
        omegaapg1, apg1 = parse_avepp(apg, max_freq)
        aac = np.append(aac, np.average(apc1))
        aas = np.append(aas, np.average(aps1))
        aag = np.append(aag, np.average(apg1))
        x = np.append(x, max_freq)
    plt.figure(figsize=(16, 16))
    plt.legend(loc="lower right")
    plt.plot(x, aac*cz**2, label="ap_alpha")
    plt.plot(x, aas*sz**2, label="ap_beta")
    plt.plot(x, aag*gz**2, label="ap_gamma")
    print x[130], aac[130]*cz**2, aas[130]*sz**2, aag[130]*gz**2
    print x[329], aac[329]*cz**2, aas[329]*sz**2, aag[329]*gz**2
    print x[100], aac[100]*cz**2, aas[100]*sz**2, aag[100]*gz**2
    print x[280], aac[280]*cz**2, aas[280]*sz**2, aag[280]*gz**2
    plt.legend(loc="upper right")
    #plt.savefig("tst_plot.pdf")

run()
plt.show()
