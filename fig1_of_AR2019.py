#!/usr/bin/env python
# This script plot the results of bin-width optimization, "reulst.txt_vec"
# to make Fig.1 for MLF annual report 2019.
# 2020 8/1 Kazuyoshi TATSUMI
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['font.family'] = 'Arial'
fig = plt.figure(figsize=(6, 40/16.0*5))
#fig.suptitle("ei42 cond09 orthotope data")
fig.suptitle("ei24 cond07 orthotope data")
#fig.suptitle("17714 cond09 orthotope opt data")


def get_data(infile, m):
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
            extraFound += 1
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
    if m > 1:
        while ForFound < m - 1:
            line = f.readline()
            if line.startswith('For'):
                ForFound += 1
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


def plotter(t_n, t_m, list_n, list_m, deltas, m, ulm):
    delta = deltas[0]*deltas[1]*deltas[2]*deltas[3]
    y_n = 1.0/(list_n[1, :]*list_n[2, :]*list_n[3, :]*list_n[4, :]*delta)
    y_m = 1.0/(list_m[1, :]*list_m[2, :]*list_m[3, :]*list_m[4, :]*delta)
    ymm = y_m.shape[0]
    ax = fig.add_subplot(5, 1, 1)
    ax.scatter(
            list_n[0, m:], y_n[m:],
            marker='s', edgecolor="black", s=50, facecolors="white")
    ax.scatter(
            list_n[0, 0:m-1], y_n[0:m-1],
            marker='s', edgecolor="black", s=50, facecolors="white")
    ax.scatter(
            list_n[0, m-1:m], y_n[m-1:m],
            marker='*', edgecolor="black", s=50, facecolors="white")
    ax.scatter(
            list_m[0, m:ymm+ulm], y_m[m:ymm+ulm],
            marker='.', edgecolor="black", s=50, facecolors="black")
    ax.scatter(
            list_m[0, 0:m-1], y_m[0:m-1],
            marker='.', edgecolor="black", s=50, facecolors="black")
    ax.set_xlim(0, list_m[0, 0]*1.02)
    ax.set_ylim(-0.03*y_m[ulm-1], y_m[ulm-1] * 1.05)
    plt.gca().ticklabel_format(style="sci", scilimits=(0, 0), axis="x")
    plt.gca().ticklabel_format(style="sci", scilimits=(0, 0), axis="y")
    ax.set_xlabel('1 / total counts')
    ax.set_ylabel('1 / delta (rlu^-3meV^-1)')
    for ip, dp in enumerate(deltas):
        y_n = list_n[1+ip, :]*dp
        y_m = list_m[1+ip, :]*dp
        ax = fig.add_subplot(5, 1, 2+ip)
        ax.scatter(
                t_n[m:], y_n[m:], label="each n",
                marker='s', edgecolor="black", s=50, facecolors="white")
        ax.scatter(
                t_n[0:m-1], y_n[0:m-1], label="each n",
                marker='s', edgecolor="black", s=50, facecolors="white")
        ax.scatter(
                t_n[m-1:m], y_n[m-1:m], label="each n",
                marker='*', edgecolor="black", s=50, facecolors="white")
        ax.scatter(
                t_m[m:ymm+ulm], y_m[m:ymm+ulm], label="each n",
                marker='.', edgecolor="black", s=50, facecolors="black")
        ax.scatter(
                t_m[0:m-1], y_m[0:m-1], label="each n",
                marker='.', edgecolor="black", s=50, facecolors="black")
        ax.set_xlim(0, max(t_m[0:ymm+ulm])*1.02)
        ax.set_ylim(0, max([np.max(y_m), np.max(y_n)])*1.05)
        ax.set_xlabel('measurement time (h)')
        ax.set_ylabel('bin width (rlu)')
        #if ip == 0:
        #    ax.set_yticks(np.arange(0, 3)*0.02)
        if ip == 3:
        #    ax.set_yticks(np.arange(0, 5)*0.2)
            ax.set_ylabel('bin width (meV)')


def get_measurement_time(xlist_n, xlist_m, hs):
    nln = len(xlist_n)
    t_n = list(range(1, (nln + 1)))
    t_m = list(range(1, (nln + 1)))
    for xm in xlist_m[nln:]:
        t_m.append((xlist_m[0] / xm))
    t_n = np.array(t_n)*hs
    t_m = np.array(t_m)*hs
    print(t_n)
    print(t_m)
    return(t_n, t_m)


def run():
    infile = "result.txt_vec"

    #dqx = 0.025  # rlu, for 17714
    #dqy = 0.025  # rlu, for 17714
    #dqz = 0.025  # rlu, for 17714
    #de  = 0.5    # meV, for 17714

    #dqx = 0.0125 # rlu, for Ei42
    #dqy = 0.025  # rlu, for Ei42
    #dqz = 0.050  # rlu, for Ei42
    #de  = 0.2    # meV, for Ei42

    dqx = 0.01   # rlu, for Ei24
    dqy = 0.01   # rlu, for Ei24
    dqz = 0.04   # rlu, for Ei24
    de  = 0.08   # meV, for Ei24


    m = 4      # sereial number of extraploation result to be plotted
    ulm = -2   # upper limit for eliminating the results of extrapolation data in the plot
    deltas = np.array([dqx, dqy, dqz, de])
    #hourperallanglescan = 2 * 146.0 / 60  # for Ei42
    #hourperallanglescan = 1.0             # for 17714
    hourperallanglescan = 4 * 11.0 / 60    # for Ei2
    list_n, list_m = get_data(infile, m)
    print(list_n.shape)
    t_n, t_m = get_measurement_time(
            list_n[0, :], list_m[0, :], hourperallanglescan
            )
    plotter(t_n, t_m, list_n, list_m, deltas, m, ulm)
    #plt.show()
    plt.savefig("ei24_cond07_orthotope_data_for_fig1.pdf")
    #plt.savefig("ei42_cond09_orthotope_data_for_fig1.pdf")
    #plt.savefig("no17714_cond09_orthotope_opt_data_for_fig1.pdf")

run()
