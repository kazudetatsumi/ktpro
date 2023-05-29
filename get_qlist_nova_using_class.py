#!/usr/bin/env python
import sys
from get_qlist_nova_class import get_qlist as gq


def run(binw, runno):
    fn = "0" + binw.split(".")[1]
    fh = "./srlz/" + fn + "io/run" + runno + "_"
    pklfile = "./srlz/" + fn + "/run" + runno + ".pkl"
    prj = gq(pklfile=pklfile)
    prj.read_pkl()
    for i in range(10, 11):
        qmin = 0.1 + 0.1*i
        qmax = qmin + 0.1
        spectrafile = fh+"_{:.2f}-{:.2f}_spectra.pkl".format(qmin, qmax)
        prj.spect(qmin, qmax, prj.dataset)
        prj.save_spectra(spectrafile)

if len(sys.argv) >= 2:
    runno = sys.argv[1]
else:
    runno = "6205"
bin = "0.000025"
run(bin, runno)
