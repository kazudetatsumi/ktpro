#!/usr/bin/env python
import sys
sys.path.append("/home/kazu/ktpro")
import qens_check_divided_spectra_class as qcdsc

def samplerun():
    if len(sys.argv) >= 2:
        runno = sys.argv[1]
    else:
        print("using default runno 6202")
        runno = "6202"
    if len(sys.argv) >= 3:
        dirname = sys.argv[2]
    else:
        print("using default dirname 0000025io")
        dirname = "0000025io"
    if len(sys.argv) >= 4:
        frac = sys.argv[3]
    else:
        print("using default frac """)
        frac = ""

    head = "/home/kazu/desktop/210108/Tatsumi/srlz/"
    histofile = head + dirname + "/run" + runno + "united_spectra.pkl"
    histifile = head + dirname + "/run" + runno + "united_monispectra.pkl"
    dividedkdefile = "./qens_kde_o_divided_by_i_" + runno + ".pkl"
    proj = qcdsc.check_divided_spectra(histofile, histifile, dividedkdefile)
    proj.plotter()


samplerun()
