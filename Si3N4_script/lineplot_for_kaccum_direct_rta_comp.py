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


def eachplot(sn, comp, omegar,  dkaccumr,  omegad, dkaccumd, max_freq, st1, st2, w1, w2):
    #plt.subplot(2, 1, sn, aspect=0.3)
    plt.subplot(2, 1, sn)
    #plt.subplot(2, 1, sn, aspect=0.45)
    plt.subplots_adjust(hspace=0, wspace=0)
    plt.plot(omegad, dkaccumd, label=comp + "_rta", linestyle='-', color='k', linewidth=w1)
    #ax[sn].plot(omegad, dkaccumd)
    #ax[1, 1].fill(omegar, dkaccumr, facecolor='lightgray', linewidth=0.1)
    #plt.plot(omegad, dkaccumd, label=comp + "_direct", linestyle='-', color='k', linewidth=w2)
    plt.fill(omegar, dkaccumr, facecolor='gray', linewidth=w2)
    #plt.xlim(0, max_freq)
    #ax[sn].set_ylim(0, 35)
    #ax[sn].set_xticks([0, 5, 10, 15, 20, 25, 30])
    #fig.subplots_adjust(hspace=0, wspace=0)
    plt.xlim(0, max_freq)
    if sn == 2:
       plt.ylim(0, 15)
       plt.yticks([0, 10])
    else:
       x = range(0, max_freq, 5)
       plt.xticks(x, " ")
       plt.ylim(0, 35)
       plt.yticks(range(0, 35, 10))
    #   ax[sn].set_yticks([0, 10, 20 ,30 ])
    #   ax[sn].set_aspect(0.35)

def run():
    sdir = "/home/kazu//bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/"
    Temp = 300
    max_freq = 33
    gs = sdir + 'noiso/kaccum_m101026_dense.dat'
    gsd = sdir + 'direct-sol/kaccum_m101026_dense.dat'
    sv = 143.78


    omegas, kaccums, dkaccums = parse_kaccum(gs)
    omegasd, kaccumsd, dkaccumsd = parse_kaccum(gsd)
    plt.figure(figsize=(4.5, 2.5))
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

    eachplot(1, "zz", omegas, dkaccums[:, 2], omegasd, dkaccumsd[:, 2], max_freq, (0, (1.7, 1.7)), (0, (6, 1, 1, 1)), 1.0, 0.0)
    eachplot(2, "xx", omegas, dkaccums[:, 0], omegasd, dkaccumsd[:, 0], max_freq, (0, ()), (0, (5, 2)), 1.0, 0.0)
    #plt.subplots_adjust(hspace=0, wspace=0)


    plt.savefig("tst_plot.pdf")


run()
plt.show()
