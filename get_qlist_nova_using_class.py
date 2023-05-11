#!/usr/bin/env python
from get_qlist_nova_class import get_qlist as gq


def run(binw, runno):
    fn = "0" + binw.split(".")[1]
    fh = "./srlz/" + fn + "io/run" + runno + "_"
    pklfile = "./srlz/" + fn + "/run" + runno + ".pkl"
    prj = gq(pklfile=pklfile)
    prj.read_pkl()
    for i in range(0, 1):
        qmin = 0.1 + 0.1*i
        qmax = qmin + 0.1
        spectrafile = fh+"_{:.1f}-{:.1f}_spectra.pkl".format(qmin, qmax)
        prj.spect(qmin, qmax, prj.dataset)
        prj.save_spectra(spectrafile)


bin = "0.000025"
runno = "6207"
run(bin, runno)
