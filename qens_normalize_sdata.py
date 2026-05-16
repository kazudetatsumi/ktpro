#!/usr/bin/env python
import numpy as np
import sys
import matplotlib.pyplot as plt
sys.path.append("/home/kazu/ktpro")
from get_qlist_nova_class import get_qlist as gq


def eachdata(pklfile):
    prj = gq(pklfile=pklfile)
    prj.read_pkl()
    return prj.dataset


def run():
    sprefix = "/home/kazu/desktop/210108/Tatsumi/srlz/0000001/"
    oprefix = "/home/kazu/desktop/210108/Tatsumi/srlz/0000001io/"
    savefix = "/home/kazu/desktop/210108/Tatsumi/winparam_exam/160_1_0000001io_Boxcar/n200000/skde/"
    runno = "6202"
    #variables = [0.5, 0.02, 0.01, 0.01]
    qmins = np.arange(0.5, 0.6, 0.1)
    qmaxs = qmins + 0.1
    for qmin, qmax in zip(qmins, qmaxs):
        sfile = sprefix + "run" + runno + "s_" + "{:.1f}".format(qmin) + "-" "{:.1f}".format(qmax) + ".pkl"
        ofile = oprefix + "run" + runno + "__" + "{:.1f}".format(qmin) + "-" "{:.1f}".format(qmax) + "_spectra.pkl"
        spectrafile =  savefix + "run" + runno + "united_spectra.pkl"
        sdataset = eachdata(sfile)
        odataset = eachdata(ofile)
        print(np.max(sdataset['spectra']))
        sy_normalized = sdataset['spectra'] / np.max(sdataset['spectra'][40000:60000]) * np.max(odataset['spectra'][40000:60000])
        sy_normalized = sy_normalized.astype(int)
        plt.plot(sdataset['energy'], sy_normalized)
        plt.plot(odataset['energy'], odataset['spectra'])
        plt.yscale('log')
        plt.show()
        prj = gq()
        prj.spectra = sy_normalized
        prj.ene = sdataset['energy']
        prj.save_spectra(spectrafile)


run()


