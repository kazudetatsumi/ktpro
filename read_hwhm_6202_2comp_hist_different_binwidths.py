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
plt.rcParams['font.size'] = 11
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
    #bins = ["000010io"]
    #bins = ["0000003io"]
    Ms=None
    pws=None
    pfs=None
    channels=None
    bins = "0000025io"
    fracs = ["", "0875", "075", "0625",  "05", "0375", "025", "0125"]
    tcount = np.array([18970., 16098., 13380., 10621., 7794., 5318., 2826. ])
    ndes = ["0.0010", "0.00025", "0.000025"]
    ylims = [[0., 0.11], [0., 0.012]]
    eachnde(Ms, pws, pfs, channels, bins, fracs, tcount, ndes, ylims)


def eachnde(Ms, pws, pfs, channels, bins, fracs, tcount, ndes, ylims):
    fig = plt.figure(figsize=(16, 8))
    for inde, nde in enumerate(ndes):
        prj = rhc.read_hwhm(Ms, pws, pfs, channels, bins, fracs,
                     prefix=home+"/desktop/210108/Tatsumi/winparam_exam/test_reconst_python",
                     numlore=2, onlyhist=True, nde=nde)
        prj.create_array_hist()
        prj.data = prj.hwhms.squeeze()
        prj.data = prj.data[1:,:]
        prj.data[:, 0] = tcount[:]/tcount[0]
        text = "~"+prj.prefix.split('Tatsumi')[1]+"/"+bins+"_*_"+nde
        plotter(fig, prj.data, inde, ylims, text=text)

    plt.subplots_adjust(hspace=0.0)
    plt.show()


#### I took the plotter method  from qens_fit_kde_hist_class and modified it.
def plotter(fig, data, inde, ylims,  text=None, c=['lightblue', 'blue'],
            loc=None, bbox_to_anchor=None):
    for idata in [0, 1]:
        ax = fig.add_subplot(2, 3, 1+inde+3*idata)
        ax.errorbar(data[:, 0], data[:, 1+2*idata],
                    yerr=data[:, 2+2*idata], marker="o", label='hist'+str(idata),
                    ms=4, elinewidth=1, lw=0, capsize=3, c=c[idata], mfc="none")
        ax.set_ylim(ylims[idata][0], ylims[idata][1])
        ax.legend()
        if idata == 1:
            ax.set_xlabel('Fraction of event data used')
            ax.tick_params(direction='in', right=True, top=True, labelbottom=True)
        else:
            ax.tick_params(direction='in', right=True, top=True, labelbottom=False)
            if text:
                ax.text(0, 0.12, text)
    ax.set_ylabel('HWHM (meV)')
    #if loc is not None and bbox_to_anchor is not None:
    #    ax.legend(loc=loc, bbox_to_anchor=bbox_to_anchor)


samplerun()
