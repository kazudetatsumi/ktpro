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
                if count > 1:
                    break
                print("%d %s" % (count, line.strip()))
                continue
            if line.strip() == "":
                continue
            if count == 1:
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
    plt.plot(omega, dos, label=phase + "_dos", linewidth=1, color="k")
    plt.ylim(0, 5.4)
    plt.yticks([0, 1, 2, 3, 4, 5])
    if sn > 1:
        plt.yticks([0, 1, 2, 3, 4, 5], " ")
    plt.xlim(0, max_freq)
    x = range(0, max_freq, 5)
    plt.xticks(x, " ")


def eachplot_kappa(sn, phase, omega, kaccum, dkaccum, numr, max_freq):
    plt.subplot(numr, 3, sn)
    # plt.title("kaccum_for_" + phase)
    # plt.plot(omega,kaccum[:,0],label=phase + "_kxx")
    # plt.plot(omega,kaccum[:,2],label=phase + "_kzz")
    plt.plot(omega, dkaccum[:, 0], label=phase + "_dkxx", color="k",  linewidth = 1.4, linestyle=(0, (1.5, 2.5)))
    if sn != 6:
        plt.plot(omega, dkaccum[:, 2], label=phase + "_dkzz", linewidth=1.0,  color="k")
    plt.ylim(0, 27)
    plt.yticks([0, 5, 10, 15, 20, 25])
    if sn > 4:
        plt.yticks([0, 5, 10, 15, 20, 25], " ")
    plt.xlim(0, max_freq)
    x = range(0, max_freq, 5)
    plt.xticks(x, " ")


def eachplot_gv(sn, phase, omega, dkaccum, numr, max_freq):
    plt.subplot(numr, 3, sn)
    plt.subplots_adjust(hspace=0, wspace=0)
    # plt.title("gv_for_" + phase)
    plt.plot(omega, dkaccum[:, 0], label=phase + "_dkxx",  color="k", linewidth = 1.4, linestyle=(0, (1.5, 2.5)))
    if sn != 9:
        plt.plot(omega, dkaccum[:, 2], label=phase + "_dkzz", color="k",  linewidth=1.0)
    plt.ylim(-3, 11)
    plt.yticks([0, 2, 4, 6, 8, 10])
    if sn > 7:
        plt.yticks([0, 2, 4, 6, 8, 10], " ")
    plt.xlim(0, max_freq)
    x = range(0, max_freq, 5)
    plt.xticks(x, " ")


def eachplot_tau(sn, phase, omega, gamma, xmin, xmax, ymin, ymax, title, numr):
    plt.subplot(numr, 3, sn)
    plt.scatter(omega, gamma, marker="o", color="k",  s=2, edgecolors="none",
                label=phase + "_" + title)
    plt.ylim(ymin, ymax)
    plt.yticks([0, 20, 40, 60, 80])
    if sn > 10:
        plt.yticks([0, 20, 40, 60, 80], " ")
    plt.xlim(xmin, xmax)
    plt.xlabel("Frequency (THz)")


def rebin(x, y): 
    rebinx = np.arange(0, 35, 0.01)
    rebiny = np.zeros((rebinx.shape[0],y.shape[1]))
    for i in range(0, rebiny.shape[0]):
        for j in range(0, x.shape[0]-1):
            if rebinx[i] >= x[j] and rebinx[i] < x[j+1]:
                rebiny[i] = y[j] + (y[j+1] - y[j]) / (x[j+1] - x[j]) * (rebinx[i] - x[j])
    return(rebinx, rebiny)


def run():
    headdir = "/home/kazu/for_figure_dos_kappa_gv_tau"
    cdir = headdir + "/asi3n4/"
    sdir = headdir + "/bsi3n4/"
    gdir = headdir + "/gsi3n4/"
    wadir = "/home/kazu/waln/phono3py_332_fc2_443_sym_monk/"
    zadir = "/home/kazu/zbaln/phonopy_444/"
    Temp = 300
    numr = 7
    max_freq = 33
    dosc = np.loadtxt(cdir + 'total_dos_m292935.dat',
                      comments='#', dtype='float')
    doss = np.loadtxt(sdir + 'total_dos_m292967.dat',
                      comments='#', dtype='float')
    dosg = np.loadtxt(gdir + 'total_dos.dat',
                      comments='#', dtype='float')
    doswa = np.loadtxt(wadir + 'total_dos_m727239.dat',
                      comments='#', dtype='float')
    dosza = np.loadtxt(zadir + 'total_dos_m727272.dat',
                      comments='#', dtype='float')
    gc = cdir + 'noiso/kaccum_m101014_dense.dat'
    gs = sdir + 'noiso/kaccum_m101026_dense.dat'
    gg = gdir + 'noiso/kaccum_m181818_dense.dat'
    ggc = cdir + 'gvaccum_m101014_dense.dat'
    ggs = sdir + 'gvaccum_m101026_dense.dat'
    ggg = gdir + 'gvaccum_m181818_dense.dat'
    c = cdir + "noiso/kappa-m101014.noiso.hdf5"
    s = sdir + "noiso/kappa-m101026.noiso.hdf5"
    g = gdir + "noiso/kappa-m181818.hdf5"
    cv = 298.78
    sv = 143.78
    gv = 472.14
    cz = 4
    sz = 2
    gz = 2

    #omegac, kaccumc, dkaccumc = parse_kaccum(gc)
    #omegas, kaccums, dkaccums = parse_kaccum(gs)
    #omegag, kaccumg, dkaccumg = parse_kaccum(gg)
    omegagc, gvaccumc, dgvaccumc = parse_gvaccum(ggc)
    omegags, gvaccums, dgvaccums = parse_gvaccum(ggs)
    #omegagg, gvaccumg, dgvaccumg = parse_gvaccum(ggg)
    #omegac1, gammac1 = parse_gamma(c, Temp, max_freq, flag=0)
    #omegas1, gammas1 = parse_gamma(s, Temp, max_freq, flag=0)
    #omegag1, gammag1 = parse_gamma(g, Temp, max_freq, flag=0)
    #romegagc, rdgvaccumc = rebin(omegagc, dgvaccumc)
    #romegags, rdgvaccums = rebin(omegags, dgvaccums)

    plt.figure(figsize=(11, 10.050))
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

    #eachplot_dos(1, "alpha", dosc[:, 0], dosc[:, 1] / cv, numr, max_freq)
    #eachplot_dos(2, "beta", doss[:, 0], doss[:, 1] / sv, numr, max_freq)
    #eachplot_dos(3, "gamma", dosg[:, 0], dosg[:, 1]/gv*4, numr, max_freq)
    #eachplot_dos(1, "alpha", dosc[:, 0], dosc[:, 1] / cz, numr, max_freq)
    #eachplot_dos(2, "beta", doss[:, 0], doss[:, 1] / sz, numr, max_freq)
    #eachplot_dos(3, "gamma", dosg[:, 0], dosg[:, 1] / gz, numr, max_freq)

    #eachplot_kappa(4, "alpha", omegac, kaccumc, dkaccumc, numr, max_freq)
    #eachplot_kappa(5, "beta", omegas, kaccums, dkaccums, numr, max_freq)
    #eachplot_kappa(6, "gamma", omegag, kaccumg, dkaccumg, numr, max_freq)

    #eachplot_gv(7, "alpha", omegagc, dgvaccumc, numr, max_freq)
    #eachplot_gv(8, "beta", omegags, dgvaccums, numr, max_freq)
    #eachplot_gv(9, "gamma", omegagg, dgvaccumg, numr, max_freq)
    #eachplot_gv(7,"beta-alpha",romegagc,rdgvaccums-rdgvaccumc, numr, max_freq)
    plt.subplot(3,1,1)
    #plt.plot(dosc[:,0], (doss[:,1] * 2 - dosc[:,1]) / (doss[:,1] * 2 + dosc[:,1]))
    ints = np.sum(doss[:, 1])*(doss[1, 0]-doss[0, 0])
    intc = np.sum(dosc[:, 1])*(dosc[1, 0]-dosc[0, 0])
    intg = np.sum(dosg[:, 1])*(dosg[1, 0]-dosg[0, 0])
    intwa = np.sum(doswa[:, 1])*(doswa[1, 0]-doswa[0, 0])
    intza = np.sum(dosza[:, 1])*(dosza[1, 0]-dosza[0, 0])
    #plt.plot(dosc[:,0], abs((doss[:, 1]/ints - dosc[:, 1]/intc)))
    plt.plot(dosc[:,0], dosc[:, 1]/intc)
    plt.plot(doss[:,0], doss[:, 1]/ints)
    #plt.plot(dosza[:,0], abs((doswa[:, 1]/intwa - dosza[:, 1]/intza)))
    plt.subplot(3,1,2)
    plt.plot(doswa[:,0], doswa[:, 1]/intwa)
    plt.plot(dosza[:,0], dosza[:, 1]/intza)
    #plt.plot(romegagc, (rdgvaccums[:,2]-rdgvaccumc[:,2]) /  (rdgvaccums[:,2]+rdgvaccumc[:,2]))
    #plt.plot(romegagc, (rdgvaccums[:,2]-rdgvaccumc[:,2]))
    #plt.plot(romegagc, rdgvaccumc[:,2])
    #plt.plot(doss[:,0], abs((doss[:, 1]/ints - dosg[:, 1]/intg)))
    plt.subplot(3,1,3)
    plt.plot(dosc[:,0], ((doss[:, 1]/ints - dosc[:, 1]/intc))**2)
    plt.plot(dosza[:,0], ((doswa[:, 1]/intwa - dosza[:, 1]/intza))**2)
    print np.sum((doss[:1700, 1]/ints - dosc[:1700, 1]/intc)**2)
    print np.sum((doswa[:1700, 1]/intwa - dosza[:1700, 1]/intza)**2)
    #plt.plot(romegagc, (rdgvaccums[:,0]-rdgvaccumc[:,0]) /  (rdgvaccums[:,0]+rdgvaccumc[:,0]))
    #plt.plot(romegagc, (rdgvaccums[:,0]-rdgvaccumc[:,0]))
    #plt.plot(romegagc, rdgvaccums[:,2])

    #eachplot_tau(10, "alpha", omegac1, 1.0 / (4 * np.pi * gammac1),
    #             0, max_freq, 0, 100, "gamma", numr)
    #eachplot_tau(11, "beta", omegas1, 1.0 / (4 * np.pi * gammas1),
    #             0, max_freq, 0, 100, "gamma", numr)
    #eachplot_tau(12, "gamma", omegag1, 1.0 / (4 * np.pi * gammag1),
    #             0, max_freq, 0, 100, "gamma", numr)

    #plt.savefig("tst_plot.pdf")


run()
plt.show()
