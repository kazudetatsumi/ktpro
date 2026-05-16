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
    Ms = None
    pws = None
    pfs = None
    channels = None
    bins = "00001"
    fracs = ["", "075", "0625",  "05", "0375", "025", "0125"]
    tcount = np.array([18970., 16098., 13380., 10621., 7794., 5318., 2826.])
    ndes = ["0.00025", "0.0001", "0.000025"]
    ylims = [[0., 0.11], [0., 0.012]]
    eachnde(Ms, pws, pfs, channels, bins, fracs, tcount, ndes, ylims)


def eachnde(Ms, pws, pfs, channels, bins, fracs, tcount, ndes, ylims):
    fig = plt.figure(figsize=(16, 8))
    for inde, nde in enumerate(ndes):
        dn = "0"+nde.split(".")[1]+"s"
        prj = rhc.read_hwhm(Ms, pws, pfs, channels, bins, fracs,
                            prefix=home+"/desktop/210108/Tatsumi/" +
                            "winparam_exam/test_shist/"+dn,
                            numlore=2, onlyhist=True, nde=nde)
        prj.create_array_shist()
        prj.data = prj.hwhms.squeeze()
        prj.data[:, 0] = tcount[:]/tcount[0]
        text = "~"+prj.prefix.split('Tatsumi')[1]+"/"+bins+"_*_"+nde
        plotter(fig, prj.data, prj.negs, inde, ylims, text=text)

    plt.subplots_adjust(hspace=0.0)
    plt.show()


#### I took the plotter method  from qens_fit_kde_hist_class and modified it.
def plotter(fig, data, negs, inde, ylims,  text=None, c=['lightblue', 'blue'],
            loc=None, bbox_to_anchor=None):
    for idata in [0, 1]:
        ax = fig.add_subplot(2, 3, 1+inde+3*idata)
        ax.errorbar(data[~negs, 0], data[~negs, 1+2*idata],
                    yerr=data[~negs, 2+2*idata], marker="o", label='hist'+str(idata),
                    ms=4, elinewidth=1, lw=0, capsize=3, c=c[idata], mfc="none")
        ax.errorbar(data[negs, 0], data[negs, 1+2*idata],
                    yerr=data[negs, 2+2*idata], marker="o", label='hist'+str(idata)+'neg',
                    ms=4, elinewidth=1, lw=0, capsize=3, c='gray', mfc="none")
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
