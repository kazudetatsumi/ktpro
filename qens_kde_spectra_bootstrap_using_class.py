#!/usr/bin/env python
# hwhms calculated with boot-strap sampled kde profiles
# Kazuyoshi Tatsumi 2022/04/07
import os
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
home=os.path.expanduser("~")
from qens_kde_results_odata_divided_by_idata_class \
 import odata_divided_by_idata as div


def run():
    elim = [-0.03, 0.07]
    odfile = home+"/desktop/210108/Tatsumi/winparam_exam/" +\
        "test_python_de_0.000025/160_1_0000025io_Boxcar_/" +\
        "qens_run6204united_kde_results_on_data_qsel.pkl"
    idfile = home+"/desktop/210108/Tatsumi/pickles/" +\
        "0000025io/qens_kde_results_on_idata_6204.pkl"
    ofile = home+"/desktop/210108/Tatsumi/winparam_exam/" +\
        "test_python_de_0.000025/160_1_0000025io_Boxcar_/" +\
        "qens_run6202united_kde_results_on_data_qsel.pkl"
    ifile = home+"/desktop/210108/Tatsumi/pickles/" +\
        "0000025io/qens_kde_results_on_idata_6202.pkl"
    fig = plt.figure(figsize=(12, 12))
    prj = div(odfile, idfile, bootstrap=True)
    prj.get_data()
    plotter(prj, fig, elim, 0)
    prj = div(ofile, ifile, bootstrap=True)
    prj.get_data()
    plotter(prj, fig, elim, 1)
    plt.show()


def plotter(prj, fig,elim, icol):
            text=prj.ofile.split("winparam_exam")[1]
            mask=np.where((prj.xo >= elim[0]) & (prj.xo <= elim[1]))
            ax = fig.add_subplot(3, 2, 1+icol)
            for i in range(0, 100):
                ax.plot(prj.xo[mask], prj.yo_ssvk[6][i][mask], lw=0.2)
            ax.set_title(text,fontsize=7.5)
            ax.set_yscale('log')
            ax = fig.add_subplot(3, 2, 3+icol)
            for i in range(0, 100):
                ax.plot(prj.xo[mask], prj.yi_ssvk_ip[i][mask], lw=0.2)
            ax.set_yscale('log')
            ax = fig.add_subplot(3, 2, 5+icol)
            for i in range(0, 100):
                ax.plot(prj.xo[mask], prj.y_ssvk[i][mask], lw=0.2)
            ax.set_yscale('log')
            plt.subplots_adjust(hspace=0.0)



run()
