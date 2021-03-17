#!/usr/bin/env python
import sys
sys.path.append("/home/kazu/ktpro")
import mergetxt as mg
import alpha_17714_Ei42_Ei24 as a1EE
import get_all_totalcount_class as gatc
import numpy as np
import matplotlib.pyplot as plt


class SMerge_Txt_Count(mg.Merge_Txt, gatc.Get_All_TCount):
    def __init__(self, yfile, workdir, tail):
        self.yfile = yfile
        self.workdir = workdir
        self.tail = tail

    def merge(self, alpha=1.0):
        super(SMerge_Txt_Count, self).get_dirnames()
        super(SMerge_Txt_Count, self).get_tcounts()
        xdirnames, xvals = np.array(self.dirnames), np.array(self.counts)*alpha
        ydirnames, yvals = super(SMerge_Txt_Count, self).read_txt(self.yfile)
        print(xdirnames)
        print(xvals)
        print(ydirnames)
        print(yvals)

        mask = np.isin(xdirnames, ydirnames)
        self.x = xvals[mask]
        xdir = xdirnames[mask]
        mask = np.isin(ydirnames, xdirnames)
        ydir = ydirnames[mask]
        self.y = yvals[mask]
        print("-----check elements ordered correctly-----")
        print(self.x)
        print(self.y)
        print(xdir)
        print(ydir)

    def plotdata(self, label, color):
        title = "max mask fractions wrt count in common space"
        super(SMerge_Txt_Count, self).create_fig(title=title)
        super(SMerge_Txt_Count, self).plotter(label, color)


def samplerun():
    workdir = "/home/kazu/desktop/200204/fine/hourbyhour/" +\
              "ortho_opt_without_mask"
    tail = "h"
    yfile = "/home/kazu/desktop/200204/fine/hourbyhour/add_random_mask/" +\
            "maxfrac.txt"
    proj = SMerge_Txt_Count(yfile, workdir, tail)
    proj.merge()
    proj.plotter(r'#1($\alpha$=0)', 'red', alpha=a1EE.alpha('17714'))

    yfile = "/home/kazu/desktop/200204/fine/hourbyhour/add_random_mask/" +\
            "maxfrac_condparam09.txt"
    proj = SMerge_Txt_Count(yfile, workdir, tail)
    proj.merge()
    proj.plotter(r'#1($\alpha$=0.9)', 'brown', alpha=a1EE.alpha('17714'))

    workdir = "/home/kazu/desktop/200522/Ei42/veryfineq/"
    tail = "m"
    yfile = "/home/kazu/desktop/200522/Ei42/veryfineq/add_random_mask/" +\
            "maxfrac.txt"
    proj = SMerge_Txt_Count(yfile, workdir, tail)
    proj.merge()
    proj.plotter(r'#2($\alpha$=0)', 'limegreen', alpha=a1EE.alpha('Ei42'))

    yfile = "/home/kazu/desktop/200522/Ei42/veryfineq/add_random_mask/" +\
            "maxfrac_condparam09.txt"
    proj = SMerge_Txt_Count(yfile, workdir, tail)
    proj.merge()
    proj.plotter(r'#2($\alpha$=0.9)', 'darkgreen', alpha=a1EE.alpha('Ei42'))

    workdir = "/home/kazu/desktop/200522/Ei24/fineq/"
    tail = "m"
    yfile = "/home/kazu/desktop/200522/Ei24/fineq/add_random_mask/" +\
            "maxfrac.txt"
    proj = SMerge_Txt_Count(yfile, workdir, tail)
    proj.merge()
    proj.plotter(r'#3($\alpha$=0)', 'blue', alpha=a1EE.alpha('Ei24'))

    yfile = "/home/kazu/desktop/200522/Ei24/fineq/add_random_mask/" +\
            "maxfrac_condparam07.txt"
    proj = SMerge_Txt_Count(yfile, workdir, tail)
    proj.merge()
    proj.plotter(r'#3($\alpha$=0.7)', 'k', alpha=a1EE.alpha('Ei24'))

    plt.xlabel('count in common space')
    plt.ylabel('upper limit for additional mask volume fraction (%)')
    plt.tick_params(top=True, right=True, direction='in', which='both')
    plt.xscale('log')

    plt.legend()

    plt.show()


samplerun()
