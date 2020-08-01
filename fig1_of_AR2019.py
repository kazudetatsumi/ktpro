#!/usr/bin/env python
# This script plot the results of bin-width optimization, "reulst.txt_vec"
# to make Fig.1 for MLF annual report 2019.
# 2020 8/1 Kazuyoshi TATSUMI
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['font.family'] = 'Arial'
fig = plt.figure(figsize=(6, 40/16.0*5))
fig.suptitle("ei42 without_cond orthotope_opt")


def get_data(infile):
    f = open(infile, 'r')
    extraFound = 0
    line = f.readline()
    xlist_n = []
    ylist_n = []
    Dqxlist_n = []
    Dqylist_n = []
    Dqzlist_n = []
    Dwlist_n = []
    while extraFound == 0:
        line = f.readline()
        if line.startswith('extrapolation'):
            extraFound = 1
        else:
            values = line.split()
            if (abs(float(values[5]) - 1.693150e-06) > 1e-12):
                xlist_n.append(float(values[5]))
                ylist_n.append(float(values[6]))
                Dqxlist_n.append(int(values[1]))
                Dqylist_n.append(int(values[2]))
                Dqzlist_n.append(int(values[3]))
                Dwlist_n.append(int(values[4]))
    line = f.readline()
    xlist_m = []
    Dqxlist_m = []
    Dqylist_m = []
    Dqzlist_m = []
    Dwlist_m = []
    ForFound = 0
    while ForFound == 0:
        line = f.readline()
        if line.startswith('For'):
            ForFound = 1
        elif not line:
            break
        else:
            values = line.split()
            if (abs(float(values[5]) - 1.69315e-06) > 1e-12):
                Dqxlist_m.append(int(values[1]))
                Dqylist_m.append(int(values[2]))
                Dqzlist_m.append(int(values[3]))
                Dwlist_m.append(int(values[4]))
                xlist_m.append(float(values[5]))
    list_m = np.array([xlist_m, Dqxlist_m, Dqylist_m, Dqzlist_m, Dwlist_m])
    list_n = np.array([xlist_n, Dqxlist_n, Dqylist_n, Dqzlist_n, Dwlist_n])
    return(list_n, list_m)


def plotter(t_n, t_m, list_n, list_m, deltas):
    delta = deltas[0]*deltas[1]*deltas[2]*deltas[3]
    y_n = 1.0/(list_n[1, :]*list_n[2, :]*list_n[3, :]*list_n[4, :]*delta)
    y_m = 1.0/(list_m[1, :]*list_m[2, :]*list_m[3, :]*list_m[4, :]*delta)
    ax = fig.add_subplot(5, 1, 1)
    ax.scatter(
            list_n[0, :], y_n,
            marker='s', edgecolor="black", s=50, facecolors="white")
    ax.scatter(
            list_m[0, 1:], y_m[1:],
            marker='.', edgecolor="black", s=50, facecolors="black")
    ax.set_xlim(list_m[0, -1]-0.5e-6, list_m[0, 0]+0.5e-6)
    ax.set_ylim(0, y_m[-1] + 0.1e4)
    plt.gca().ticklabel_format(style="sci", scilimits=(0, 0), axis="x")
    plt.gca().ticklabel_format(style="sci", scilimits=(0, 0), axis="y")
    ax.set_xlabel('1 / total counts')
    ax.set_ylabel('1 / delta (rlu^-3meV^-1)')
    for ip, dp in enumerate(deltas):
        ax = fig.add_subplot(5, 1, 2+ip)
        ax.scatter(
                t_n[1:], list_n[1+ip, 1:]*dp, label="each n",
                marker='s', edgecolor="black", s=50, facecolors="white")
        ax.scatter(
                t_m[1:], list_m[1+ip, 1:]*dp, label="each n",
                marker='.', edgecolor="black", s=50, facecolors="black")
        ax.set_ylim(0, (max(list_m[1+ip, 1:])+0.5)*dp)
        ax.set_xlabel('measurement time (h)')
        ax.set_ylabel('bin width (rlu)')
        if ip == 0:
            ax.set_yticks(np.arange(0, 3)*0.02)
        if ip == 3:
            ax.set_yticks(np.arange(0, 5)*0.2)
            ax.set_ylabel('bin width (meV)')


def get_measurement_time(xlist_n, xlist_m, hs):
    t_n = range(2, (len(xlist_n)+1)*2, 2)
    t_m = range(2, (len(xlist_n)+1)*2, 2)
    for xm in xlist_m[6:]:
        t_m.append(2 * (xlist_m[0] / xm))
    t_n = np.array(t_n)*hs
    t_m = np.array(t_m)*hs
    return(t_n, t_m)


def run():
    infile = "result.txt_vec"
    dqx = 0.0125  # rlu
    dqy = 0.025  # rlu
    dqz = 0.050  # rlu
    de = 0.2   # meV
    deltas = np.array([dqx, dqy, dqz, de])
    hourperanglescan = 146.0 / 60
    list_n, list_m = get_data(infile)
    print(list_n.shape)
    t_n, t_m = get_measurement_time(
            list_n[0, :], list_m[0, :], hourperanglescan
            )
    plotter(t_n, t_m, list_n, list_m, deltas)
    #plt.show()
    plt.savefig("ei42_condparam09_for_fig1.pdf")

run()
