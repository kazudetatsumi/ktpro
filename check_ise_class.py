#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/kazu/ktpro")
import gather_optbinindx_class as goc
import re


class SCheckiselog(goc.gather_optbinidx):
    def __init__(self, log, reflog):
        self.log = log
        self.reflog = reflog
        self.refbinidx = super(SCheckiselog, self).getbinidx(self.reflog)

    def getiselogdata(self):
        ise = []
        dcount = []
        bw = []
        for line in open(self.log):
            if not re.compile('^[A-z[]').search(line):
                values = line.split()
                if int(values[3]) <= self.refbinidx[3]*2 and\
                   int(values[4]) <= self.refbinidx[2]*2-1 and\
                   int(values[5]) <= self.refbinidx[1]*2 and\
                   int(values[6]) <= self.refbinidx[0]*2:
                    ise.append(float(values[2]))
                    dcount.append(int(float(values[8])))
                    bw.append([int(values[6]), int(values[5]), int(values[4]),
                               int(values[3])])
                if int(values[3]) == self.refbinidx[3] and\
                   int(values[4]) == self.refbinidx[2] and\
                   int(values[5]) == self.refbinidx[1] and\
                   int(values[6]) == self.refbinidx[0]:
                    self.refdcount = int(float(values[8]))
                    self.refise = float(values[2])

        self.ise = np.array(ise)
        self.dcount = np.array(dcount)
        self.bw = np.array(bw)

    def plot_ise_dcount(self):
        plt.xscale('log')
        #plt.yscale('log')
        plt.scatter(self.dcount, self.ise, marker='x')
        plt.scatter(self.refdcount, self.refise, marker='x', c='gray')
        plt.scatter(self.dcount[self.ise == np.min(self.ise)], np.min(self.ise), marker='x', c='k')


def samplerun():
    log = "./test.log"
    log = "./std-ise-6_18m.log"
    reflog = "./std-18m.log"
    prj = SCheckiselog(log, reflog)
    prj.getiselogdata()
    prj.plot_ise_dcount()

    plt.show()

#samplerun()
