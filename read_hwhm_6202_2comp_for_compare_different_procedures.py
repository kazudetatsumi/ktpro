#!/usr/bin/env python
import sys
import os
home=os.path.expanduser("~")
sys.path.append(home+"/ktpro")
import read_hwhm_class as rhc
import numpy as np
import matplotlib.pyplot as plt

params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'Arial'


def samplerun():
    if len(sys.argv) >= 2:
        if sys.argv[1] == "fixbg":
            fixbg = True
        else:
            fixbg = False

    if len(sys.argv) >= 4:
        ymin = float(sys.argv[2])
        ymax = float(sys.argv[3])
    else:
        ymin = None
        ymax = None

    #Ms = [80, 160, 320, 640, 1280, 2560]
    #pws = [0.0625, 0.125, 0.25, 0.5, 1, 2, 5]
    #pfs = ["Boxcar", "Gauss"]
    #channels = ["", "channel"]
    #bins = ["000025io", "000010io", "0000025io", "0000003io"]
    #fracs = ["", "0875", "075", "05", "0375"]
    Ms = [160]
    pws = [0.25]
    pfs = ["Boxcar"]
    channels = [""]
    #bins = ["000010io"]
    #bins = ["0000003io"]
    bins = ["0000025io"]
    fracs = ["", "0875", "075", "0625",  "05", "0375", "025", "0125"]
    tcount = np.array([18970., 16098., 13380., 10621., 7794., 5318., 2826. ])
    #fracs = ["", "0875", "075", "0625",  "05", "0375"]
    prj1 = rhc.read_hwhm(Ms, pws, pfs, channels, bins, fracs,
                         #prefix=home+"/desktop/210108/Tatsumi/winparam_exam/test_python_de8",
                         prefix=home+"/desktop/210108/Tatsumi/winparam_exam/",
                         numlore=2, fixbg=False)
    prj1.create_array()
    prj1.data = prj1.hwhms.squeeze()
    prj1.data = prj1.data[1:,:]
    prj1.data[:, 0] = tcount[:]/tcount[0]

    #bins = ["000025io"]
    bins = ["0000025io"]
    pws = [0.25]
    prj2 = rhc.read_hwhm(Ms, pws, pfs, channels, bins, fracs,
                         prefix=home+"/desktop/210108/Tatsumi/winparam_exam/direct_de_idata/",
                         #prefix=home+"/desktop/210108/Tatsumi/winparam_exam/test_python_de10",
                         numlore=2, fixbg=False)
    prj2.create_array()
    prj2.data = prj2.hwhms.squeeze()
    prj2.data = prj2.data[1:,:]
    prj2.data[:, 0] = tcount[:]/tcount[0]

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(2, 2, 1)
    text = "~/"+prj1.prefix.split('Tatsumi')[1]+"@"+prj1.bins[0]+"_"+str(prj1.pws[0])
    plotter(ax, prj1.data, ymin=-0.002, ymax=0.11, text=text)
    ax = fig.add_subplot(2, 2, 3)
    plotter(ax, prj1.data, isbottom=True, ymin=-0.001, ymax=0.0125)

    ax = fig.add_subplot(2, 2, 2)
    text = "~/"+prj2.prefix.split('Tatsumi')[1]+"@"+prj2.bins[0]+"_"+str(prj2.pws[0])
    plotter(ax, prj2.data, ymin=-0.002, ymax=0.11, text=text)
    ax = fig.add_subplot(2, 2, 4)
    plotter(ax, prj2.data, isbottom=True, ymin=-0.001, ymax=0.0125)

    #print(prj.hwhms.squeeze())
    #print(prj.hwhms.squeeze().shape)
    plt.subplots_adjust(hspace=0.0)
    #plt.savefig("hwhm_for_qbsf2021.pdf")
    plt.show()


#### I took the plotter method  from qens_fit_kde_hist_class and modified it.
def plotter(ax, data, dataname=None, isbottom=False, ymin=None, ymax=None, text=None,
            c=['pink', 'red', 'lightblue', 'blue'],
            #c=['gray', 'k', 'gray', 'k'],
            loc=None, bbox_to_anchor=None):
    ax.errorbar(data[:, 0], data[:, 1],
                yerr=data[:, 2], marker="o", label='kde1',
                ms=4, elinewidth=1, lw=0, capsize=3, c=c[0], mfc="none")
    ax.errorbar(data[:, 0], data[:, 3],
                yerr=data[:, 4], marker="o", label='kde2',
                ms=4, elinewidth=1, lw=0, capsize=3, c=c[1], mfc="none")
    ax.errorbar(data[:, 0], data[:, 5],
                yerr=data[:, 6], marker="x", label='hist1',
                ms=4, elinewidth=1, lw=0, capsize=3, c=c[2])
    ax.errorbar(data[:, 0], data[:, 7],
                yerr=data[:, 8], marker="x", label='hist2',
                ms=4, elinewidth=1, lw=0, capsize=3, c=c[3])
    if text:
        ax.text(0, 0.12, text)
    if ymin is not None and ymax is not None:
        ax.set_ylim([ymin, ymax])
    if isbottom:
        ax.set_xlabel('Fraction of event data used')
        ax.tick_params(direction='in', right=True, top=True, labelbottom=True)
        ax.set_yticks([0, 0.005, 0.01])
    else:
        ax.legend()
        ax.tick_params(direction='in', right=True, top=True, labelbottom=False)
    ax.set_ylabel('HWHM (meV)')
    #if loc is not None and bbox_to_anchor is not None:
    #    ax.legend(loc=loc, bbox_to_anchor=bbox_to_anchor)


samplerun()
