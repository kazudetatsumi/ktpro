#!/usr/bin/env python
# This script plot the results of bin-width optimization, "reulst.txt_vec"
# to make Fig.1 for MLF annual report 2019.
# 2020 8/1 Kazuyoshi TATSUMI
import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib.ticker import LogFormatterSciNotation, ScalarFormatter
params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'Arial'


class Plot_4w_Logn:
    def __init__(self, infile, deltas, m, ulm=0, num_clms=1, mnj=1,
                 islog=False, isplain=True, msize=80):
        self.infile = infile
        self.deltas = deltas
        self.m = m
        self.ulm = ulm
        self.num_clms = num_clms
        self.mnj = mnj
        self.islog = islog
        self.isplain = isplain
        self.msize = msize

    def create_fig(self):
        self.fig = plt.figure(figsize=(10, 10.5))
        self.fig.suptitle("optimal bin widths")

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
                elif not line:
                    print('No line left during skipping the extrapolation' +
                          ' results in the input file! Check m value.')
                    sys.exit()
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

    def plotter(self, clm=1, ylabel=True, alpha=1.0):
        self.get_data()
        ymm = self.list_m.shape[1]
        wlist = ["$q_x$", "$q_y$", "$q_z$", "$\omega$"]
        for ip, dp in enumerate(self.deltas):
            y_n = self.list_n[1+ip, :]*dp
            x_n = 1.0/self.list_n[0, :]*alpha
            if self.islog:
                x_m = 1.0/self.list_m[0, :]*alpha
                y_m = self.list_m[1+ip, :]*dp
                deleteidx = []
                for idx in range(1, x_n.shape[0]):
                    if x_n[idx]/x_n[idx-1] < 1.2:
                        deleteidx.append(idx)
                x_n = np.delete(x_n, deleteidx)
                x_m = np.delete(x_m, deleteidx)
                y_m = np.delete(y_m, deleteidx)
                y_n = np.delete(y_n, deleteidx)
            else:
                x_m = 1.0/self.list_m[0, :-5]*alpha
                y_m = self.list_m[1+ip, :-5]*dp
                x_m = 1.0/self.list_m[0]*alpha
                y_m = self.list_m[1+ip]*dp
            ax = self.fig.add_subplot(4, self.num_clms, clm+self.num_clms*ip)
            if self.isplain:
                ax.scatter(x_n, y_n)
                ax.scatter(x_m, y_m, marker='x')
            else:
                ax.scatter(x_n[self.m:], y_n[self.m:], marker='s',
                           edgecolor="black", s=self.msize, facecolors="white")
                ax.scatter(x_n[0:self.m-1], y_n[0:self.m-1], marker='s',
                           edgecolor="black", s=self.msize, facecolors="white")
                ax.scatter(x_n[self.m-1:self.m], y_n[self.m-1:self.m],
                           marker='*', edgecolor="black", s=self.msize,
                           facecolors="white")
                ax.scatter(x_m[self.m:ymm+self.ulm], y_m[self.m:ymm+self.ulm],
                           marker='.', edgecolor="black", s=self.msize,
                           facecolors="black")
                ax.scatter(x_m[0:self.m-1], y_m[0:self.m-1],
                           marker='.', edgecolor="black", s=self.msize,
                           facecolors="black")
                ax.plot(x_n, y_n, marker='None', color="gray",
                        linestyle="dashed", lw=2)
                ax.plot(x_m[0:ymm+self.ulm], y_m[0:ymm+self.ulm],
                        marker='None', color="k", linestyle="dotted")
            ax.tick_params(labelbottom=False)
            if isinstance(self.mnj, list):
                ax.set_ylim(0, max([np.max(y_n), np.max(y_m),
                                   dp*(self.mnj[ip]-1)]) + dp)
            else:
                ax.set_ylim(0, max([np.max(y_n), np.max(y_m),
                                    dp*(self.mnj-1)]) + dp)
            if self.islog:
                ax.set_xscale('log')
                ax.text(np.max(x_m)*10**(-0.18), max(np.max(y_n), np.max(y_m)),
                        wlist[ip], size=15)
                ax.set_xlim(np.min(x_m)*10.0**(-0.06), np.max(x_m)*10.0**0.06)
            else:
                ax.text(np.max(x_m)*0.95, max(np.max(y_n), np.max(y_m)),
                        wlist[ip], size=15)
            ax.tick_params(top=True, right=True, direction='in', which='both')
            if ylabel and ip == 1:
                ax.set_ylabel('bin width (rlu)')
            if ip == 3:
                if ylabel:
                    ax.set_ylabel('bin width (meV)')
                ax.tick_params(labelbottom=True)
            ax.set_xlabel('count')
            if isinstance(self.mnj, list):
                ax.yaxis.set_major_locator(MultipleLocator(dp*self.mnj[ip]))
            else:
                ax.yaxis.set_major_locator(MultipleLocator(dp*self.mnj))
            ax.yaxis.set_minor_locator(MultipleLocator(dp))
            ax.tick_params(length=6, which='major')
            ax.tick_params(length=3, which='minor')
            plt.subplots_adjust(wspace=0.2, hspace=0.0)
