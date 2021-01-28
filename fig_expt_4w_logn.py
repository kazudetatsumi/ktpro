#!/usr/bin/env python
# This script plot the results of bin-width optimization, "reulst.txt_vec"
# to make Fig.1 for MLF annual report 2019.
# 2020 8/1 Kazuyoshi TATSUMI
import matplotlib.pyplot as plt
import numpy as np
params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)
plt.rcParams['font.size'] = 12

plt.rcParams['font.family'] = 'Arial'
fig = plt.figure(figsize=(12, 12.5))
#fig.suptitle("ei42 cond09 orthotope data")
#fig.suptitle("ei24 cond07 orthotope data")
fig.suptitle("4w_logn_with_using_mask_info")


class Plot_4w_Logn:

    def __init__(self, infile, deltas, m, ulm, num_clms):
        self.infile = infile
        self.deltas = deltas
        self.m = m
        self.ulm = ulm
        self.num_clms = num_clms

    def get_data(self):
        f = open(self.infile, 'r')
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
        if self.m > 1:
            while ForFound < self.m - 1:
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
        self.list_m = np.array([xlist_m, Dqxlist_m, Dqylist_m, Dqzlist_m,
                                Dwlist_m])
        self.list_n = np.array([xlist_n, Dqxlist_n, Dqylist_n, Dqzlist_n,
                                Dwlist_n])

    def plotter(self, clm):
        self.get_data()
        ymm = self.list_m.shape[1]
        wlist = ["$q_x$", "$q_y$", "$q_z$", "$\omega$"]
        for ip, dp in enumerate(self.deltas):
            y_n = self.list_n[1+ip, :]*dp
            y_m = self.list_m[1+ip, :]*dp
            x_n = np.log10(1.0/self.list_n[0, :])
            x_m = np.log10(1.0/self.list_m[0, :])
            ax = fig.add_subplot(4, self.num_clms, clm+self.num_clms*ip)
            ax.scatter(x_n[self.m:], y_n[self.m:], marker='s',
                       edgecolor="black", s=50, facecolors="white")
            ax.scatter(x_n[0:self.m-1], y_n[0:self.m-1], marker='s',
                       edgecolor="black", s=50, facecolors="white")
            ax.scatter(x_m[self.m:ymm+self.ulm], y_m[self.m:ymm+self.ulm],
                       marker='.', edgecolor="black", s=50, facecolors="black")
            ax.scatter(x_m[0:self.m-1], y_m[0:self.m-1],
                       marker='.', edgecolor="black", s=50, facecolors="black")
            ax.plot(x_n, y_n, marker='None', color="gray", linestyle="dashed",
                    lw=3)
            ax.plot(x_m[0:ymm+self.ulm], y_m[0:ymm+self.ulm], marker='None',
                    color="k", linestyle="dotted")
            ax.scatter(x_n[self.m-1:self.m], y_n[self.m-1:self.m], marker='*',
                       edgecolor="black", s=50, facecolors="white")
            ax.tick_params(labelbottom=False)
            ax.text(np.max(np.log10(1.0/self.list_m[0, :ymm+self.ulm]))*0.98,
                    (max(np.max(y_n), np.max(y_m))//(dp*2)+1)*(dp*2)*0.75,
                    wlist[ip], size=15)
            ax.set_ylim(0, (max(np.max(y_n), np.max(y_m))//(dp*2)+1)*(dp*2))
            xmargin = (np.max(np.log10(1.0/self.list_m[0, :])) -
                       np.min(np.log10(1.0/self.list_m[0, :])))*0.01
            ax.set_xlim(np.min(np.log10(1.0/self.list_m[0, :]))-xmargin,
                        np.max(np.log10(1.0/self.list_m[0, :]))+xmargin)
            ax.set_yticks(np.arange(0,
                          (max(np.max(y_n), np.max(y_m))//(dp*2)+2)*(dp*2),
                          (dp*2)))
            ax.tick_params(direction='in', top=True, right=True)
            ax.set_ylabel('bin width (rlu)')
            if ip == 3:
                ax.set_ylabel('bin width (meV)')
                ax.tick_params(labelbottom=True)
                ax.set_xlabel('log10(total count)')
            plt.subplots_adjust(wspace=0.4, hspace=0.0)


def samplerun():
    num_clms = 3
    m = 3
    ulm = 0

    infile = "/home/kazu/desktop/200204/fine/hourbyhour/" +\
             "ortho_opt_without_mask/condparam09/result.txt_vec"
    deltas = np.array([0.025, 0.025, 0.025, 0.5])
    clm = 1
    projectset_17714 = Plot_4w_Logn(infile, deltas, m, ulm, num_clms)
    projectset_17714.plotter(clm)

    infile = "/home/kazu/desktop/200522/Ei42/veryfineq/" +\
             "condparam_09/result.txt_vec"
    deltas = np.array([0.0125, 0.025, 0.050, 0.2])
    clm = 2
    projectset_Ei42 = Plot_4w_Logn(infile, deltas, m, ulm, num_clms)
    projectset_Ei42.plotter(clm)

    infile = "/home/kazu/desktop/200522/Ei24/fineq/" +\
             "/condparam07/result.txt_vec"
    deltas = np.array([0.01, 0.01, 0.04, 0.08])
    clm = 3
    projectset_Ei24 = Plot_4w_Logn(infile, deltas, m, ulm, num_clms)
    projectset_Ei24.plotter(clm)


samplerun()
plt.show()
