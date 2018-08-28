#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import h5py

plt.style.use('classic')
plt.rcParams['font.family'] = 'Times New Roman'


def parse_kaccum(filename):
    data = []
    with open(filename) as f:
        count = 0
        for line in f:
            if line[0] == '#':
                count += 1
                if count > 21:
                    break
                print("%d %s" % (count, line.strip()))
                continue
            if line.strip() == "":
                continue
            if count == 21:
                data.append([float(x) for x in line.split()])
    kaccum = np.array(data)[:, 1:4]
    dkaccum = np.array(data)[:, 7:10]
    omega = np.array(data)[:, 0]
    return (omega, kaccum, dkaccum)


def parse_gvaccum(filename):
    data = []
    with open(filename) as f:
        for line in f:
            if line[0] == '#':
                continue
            if line.strip() == "":
                continue
            data.append([float(x) for x in line.split()])
    gvaccum = np.array(data)[:, 1:4]
    dgvaccum = np.array(data)[:, 7:10]
    omega = np.array(data)[:, 0]
    return (omega, gvaccum, dgvaccum)


def parse_gamma(filename, temp, max_freq, flag=0):
    freqs = []
    mode_prop = []
    f = h5py.File(filename, 'r')
    temperature = f["temperature"].value
    i = 0
    for t in temperature:
        if abs(t - temp) < 1e-8:
            tindex = i
            break
        i += 1
    print("parse_gamma T=%f" % f['temperature'][tindex])
    gamma = f["gamma"][tindex, ]
    omega = f["frequency"][:, :]

    for freq, g in zip(omega, gamma):
        condition = freq < max_freq
        _freq = np.extract(condition, freq)
        _g = np.extract(condition, g)
        freqs += list(_freq)
        mode_prop += list(_g)
    gamma1 = np.array(mode_prop)
    omega1 = np.array(freqs)

    return omega1, gamma1


def eachplot_dos(sn, phase, omega, dos, numr, max_freq):
    plt.subplot(numr, 3, sn)
    plt.title(phase)
    plt.plot(omega, dos, label=phase + "_dos")
    plt.ylim(0, 0.075)
    plt.yticks([0, 0.025, 0.05, 0.075])
    if sn > 1:
        plt.yticks([0, 0.025, 0.05, 0.075], " ")
    plt.xlim(0, max_freq)
    x = range(0, max_freq, 5)
    plt.xticks(x, " ")


def eachplot_kappa(sn, phase, omega, kaccum, dkaccum, numr, max_freq):
    plt.subplot(numr, 3, sn)
    # plt.title("kaccum_for_" + phase)
    # plt.plot(omega,kaccum[:,0],label=phase + "_kxx")
    # plt.plot(omega,kaccum[:,2],label=phase + "_kzz")
    plt.plot(omega, dkaccum[:, 0] * 10, label=phase + "_dkxx")
    if sn != 6:
        plt.plot(omega, dkaccum[:, 2] * 10, label=phase + "_dkzz")
    plt.ylim(0, 250)
    plt.yticks([0, 100, 200])
    if sn > 4:
        plt.yticks([0, 100, 200], " ")
    plt.xlim(0, max_freq)
    x = range(0, max_freq, 5)
    plt.xticks(x, " ")


def eachplot_gv(sn, phase, omega, dkaccum, numr, max_freq):
    plt.subplot(numr, 3, sn)
    plt.subplots_adjust(hspace=0, wspace=0)
    # plt.title("gv_for_" + phase)
    plt.plot(omega, dkaccum[:, 0], label=phase + "_dkxx")
    if sn != 9:
        plt.plot(omega, dkaccum[:, 2], label=phase + "_dkzz")
    plt.ylim(0, 10)
    plt.yticks([0, 5, 10])
    if sn > 7:
        plt.yticks([0, 5, 10], " ")
    plt.xlim(0, max_freq)
    x = range(0, max_freq, 5)
    plt.xticks(x, " ")


def eachplot_tau(sn, phase, omega, gamma, xmin, xmax, ymin, ymax, title, numr):
    plt.subplot(numr, 3, sn)
    plt.scatter(omega, gamma, marker="o", color="k",  s=2, edgecolors="none",
                label=phase + "_" + title)
    plt.ylim(ymin, ymax)
    plt.yticks([0, 100, 200])
    if sn > 10:
        plt.yticks([0, 100, 200], " ")
    plt.xlim(xmin, xmax)
    plt.xlabel("Frequency (THz)")


def run():
    cdir = "./asi3n4/"
    sdir = "./bsi3n4/"
    gdir = "./gsi3n4/"
    Temp = 300
    numr = 7
    max_freq = 33
    dosc = np.loadtxt(cdir + 'total_dos_m292935.dat',
                      comments='#', dtype='float')
    doss = np.loadtxt(sdir + 'total_dos_m292967.dat',
                      comments='#', dtype='float')
    dosg = np.loadtxt(gdir + 'total_dos.dat',
                      comments='#', dtype='float')
    gc = cdir + 'noiso/kaccum_m101014.dat'
    gs = sdir + 'noiso/kaccum_m101026.dat'
    gg = gdir + 'noiso/kaccum_m181818.dat'
    ggc = cdir + 'gvaccum_m101014.dat'
    ggs = sdir + 'gvaccum_m101026.dat'
    ggg = gdir + 'gvaccum_m181818.dat'
    c = cdir + "noiso/kappa-m101014.noiso.hdf5"
    s = sdir + "noiso/kappa-m101026.noiso.hdf5"
    g = gdir + "noiso/kappa-m181818.hdf5"
    cv = 298.78
    sv = 143.78
    gv = 472.14

    omegac, kaccumc, dkaccumc = parse_kaccum(gc)
    omegas, kaccums, dkaccums = parse_kaccum(gs)
    omegag, kaccumg, dkaccumg = parse_kaccum(gg)
    omegagc, gvaccumc, dgvaccumc = parse_gvaccum(ggc)
    omegags, gvaccums, dgvaccums = parse_gvaccum(ggs)
    omegagg, gvaccumg, dgvaccumg = parse_gvaccum(ggg)
    omegac1, gammac1 = parse_gamma(c, Temp, max_freq, flag=0)
    omegas1, gammas1 = parse_gamma(s, Temp, max_freq, flag=0)
    omegag1, gammag1 = parse_gamma(g, Temp, max_freq, flag=0)

    plt.figure(figsize=(11, 10.050))
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

    eachplot_dos(1, "alpha", dosc[:, 0], dosc[:, 1] / cv, numr, max_freq)
    eachplot_dos(2, "beta", doss[:, 0], doss[:, 1] / sv, numr, max_freq)
    eachplot_dos(3, "gamma", dosg[:, 0], dosg[:, 1]/gv*4, numr, max_freq)

    eachplot_kappa(4, "alpha", omegac, kaccumc, dkaccumc, numr, max_freq)
    eachplot_kappa(5, "beta", omegas, kaccums, dkaccums, numr, max_freq)
    eachplot_kappa(6, "gamma", omegag, kaccumg, dkaccumg, numr, max_freq)

    eachplot_gv(7, "alpha", omegagc, dgvaccumc, numr, max_freq)
    eachplot_gv(8, "beta", omegags, dgvaccums, numr, max_freq)
    eachplot_gv(9, "gamma", omegagg, dgvaccumg, numr, max_freq)

    eachplot_tau(10, "alpha", omegac1, 1.0 / (4 * np.pi * gammac1),
                 0, max_freq, 0, 200, "gamma", numr)
    eachplot_tau(11, "beta", omegas1, 1.0 / (4 * np.pi * gammas1),
                 0, max_freq, 0, 200, "gamma", numr)
    eachplot_tau(12, "gamma", omegag1, 1.0 / (4 * np.pi * gammag1),
                 0, max_freq, 0, 200, "gamma", numr)

    plt.savefig("tst_plot.pdf")


run()
plt.show()
