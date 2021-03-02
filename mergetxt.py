#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

class Merge_Txt:
    def __init__(self, xfile, yfile):
        self.xfile = xfile
        self.yfile = yfile

    def read_txt(self, txtfile):
        dirnames = []
        vals = []
        for line in open(txtfile):
            values = line.split()
            dirnames.append(values[0])
            vals.append(float(values[1]))
        return np.array(dirnames), np.array(vals)

    def merge(self):
        xdirnames, xvals = self.read_txt(self.xfile)
        ydirnames, yvals = self.read_txt(self.yfile)
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

    def create_fig(self, title=None):
        fig = plt.figure(figsize=(8, 8)) 
        if title is None:
            fig.suptitle("max mask fractions with respect to numbers" +
                         " of bins with the optimized bin-widths")
        else:
            fig.suptitle(title)

    def plotter(self, label, color):
        plt.scatter(self.x, self.y, label=label, color=color)




def samplerun():
    yfile = "/home/kazu/desktop/200522/Ei24/fineq/add_random_mask/maxfrac.txt"
    xfile = "/home/kazu/desktop/200522/Ei24/fineq/nbins.txt"
    proj = Merge_Txt(xfile, yfile)
    proj.merge()
    proj.create_fig()
    proj.plotter('Ei24', 'blue')

    yfile = "/home/kazu/desktop/200522/Ei24/fineq/add_random_mask/maxfrac_condparam07.txt"
    xfile = "/home/kazu/desktop/200522/Ei24/fineq/condparam07/nbins.txt"
    proj = Merge_Txt(xfile, yfile)
    proj.merge()
    proj.plotter('Ei24_cond', 'k')

    yfile = "/home/kazu/desktop/200522/Ei42/veryfineq/add_random_mask/maxfrac.txt"
    xfile = "/home/kazu/desktop/200522/Ei42/veryfineq/nbins.txt"
    proj = Merge_Txt(xfile, yfile)
    proj.merge()
    proj.plotter('Ei42', 'limegreen')

    yfile = "/home/kazu/desktop/200522/Ei42/veryfineq/add_random_mask/maxfrac_condparam09.txt"
    xfile = "/home/kazu/desktop/200522/Ei42/veryfineq/condparam_09/nbins.txt"
    proj = Merge_Txt(xfile, yfile)
    proj.merge()
    proj.plotter('Ei42_cond', 'darkgreen')

    yfile = "/home/kazu/desktop/200204/fine/hourbyhour/add_random_mask/maxfrac.txt"
    xfile = "/home/kazu/desktop/200204/fine/hourbyhour/ortho_opt_without_mask/nbins.txt"
    proj = Merge_Txt(xfile, yfile)
    proj.merge()
    proj.plotter('13714', 'red')

    yfile = "/home/kazu/desktop/200204/fine/hourbyhour/add_random_mask/maxfrac_condparam09.txt"
    xfile = "/home/kazu/desktop/200204/fine/hourbyhour/ortho_opt_without_mask/condparam09/nbins.txt"
    proj = Merge_Txt(xfile, yfile)
    proj.merge()
    proj.plotter('13714_cond', 'brown')
    plt.xlabel('number of bins')
    plt.ylabel('upper limit for additional mask volume fraction')

    plt.legend()

    plt.show()


samplerun()
