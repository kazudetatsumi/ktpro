#!/usr/bin/env python
# hwhms calculated with boot-strap sampled kde profiles
# Kazuyoshi Tatsumi 2022/04/07
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/kazu/ktpro")
home=os.path.expanduser("~")
from qens_kde_hwhm_bootstrap_class import bootstrap as boot


def run():
    fracs = ["", "075", "0625", "05", "0375", "025", "0125"]
    elim = [-0.03, 0.07]
    hwhms = np.zeros((2, len(fracs)))
    shwhms = np.zeros((2, len(fracs)))
    tcount = np.array([18970., 16098., 13380., 10621., 7794., 5318., 2826.])
    x = tcount/tcount[0]
    ylims = [[0., 0.11], [0., 0.012]]
    savefile = "./results_bootstrap.pkl"
    if not os.path.exists(savefile):
        for fidx, frac in enumerate(fracs):
            print(frac)
            odfile = home+"/desktop/210108/Tatsumi/winparam_exam/" +\
                "test_python_de_0.000025/160_1_0000025io_Boxcar_/" +\
                "qens_run6204united_kde_results_on_data_qsel.pkl"
            idfile = home+"/desktop/210108/Tatsumi/pickles/" +\
                "0000025io/qens_kde_results_on_idata_6204.pkl"
            ofile = home+"/desktop/210108/Tatsumi/winparam_exam/" +\
                "test_python_de_0.000025/160_1_0000025io_Boxcar_"+frac+"/" +\
                "qens_run6202"+frac+"united_kde_results_on_data_qsel.pkl"
            ifile = home+"/desktop/210108/Tatsumi/pickles/" +\
                "0000025io/qens_kde_results_on_idata_6202"+frac+".pkl"
            prj = boot(odfile, idfile, ofile, ifile, elim)
            prj.quiet = True
            prj.divide()
            prj.preprocess()
            hwhms[:, fidx], shwhms[:, fidx] = prj.allfit()
        dataset = {}
        dataset['odfile'] = odfile
        dataset['idfile'] = idfile
        dataset['ofile'] = ofile
        dataset['ifile'] = ifile
        dataset['hwhms'] = hwhms
        dataset['shwhms'] = shwhms
        dataset['x'] = x
        with open("./results_bootstrap.pkl", 'wb') as f:
            #pickle.dump(dataset, f, -1)
            pickle.dump(dataset, f, 4)
    else:
        with open(savefile, 'rb') as f:
            #print('reading ', savefile)
            dataset = pickle.load(f)
            hwhms = dataset['hwhms']
            shwhms = dataset['shwhms']
            x = dataset['x']
    print(dataset['odfile'])
    print(dataset['idfile'])
    print(dataset['ofile'])
    print(dataset['ifile'])
    fig = plt.figure(figsize=(16, 8))
    ax = fig.add_subplot(2, 3, 1)
    ax.errorbar(x, hwhms[0, :], yerr=shwhms[0, :], marker="o", label='kde1',
                ms=4, elinewidth=1, lw=0, capsize=3, c='lightblue', mfc="none")
    ax.set_ylim(ylims[0][0], ylims[0][1])

    ax = fig.add_subplot(2, 3, 4)
    ax.errorbar(x, hwhms[1, :], yerr=shwhms[1, :], marker="o", label='kde1',
                ms=4, elinewidth=1, lw=0, capsize=3, c='blue', mfc="none")
    ax.set_ylim(ylims[1][0], ylims[1][1])
    plt.subplots_adjust(hspace=0.0)
    plt.show()


run()
