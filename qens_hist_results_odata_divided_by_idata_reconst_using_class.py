#!/usr/bin/env python
import sys
sys.path.append("/home/kazu/ktpro")
import qens_kde_results_odata_divided_by_idata_class as qkrodic

# example command:  qens_hist_results_odata_divided_by_idata_using_class.py 6204 0000025io 0.00025 0875 

def run():
    if len(sys.argv) >= 2:
        runno = sys.argv[1]
    else:
        print("using default runno 6204")
        runno = "6204"
    if len(sys.argv) >= 3:
        dirname = sys.argv[2]
    else:
        print("using default dirname 000025io")
        dirname = "000025io"
    if len(sys.argv) >= 4:
        denew = sys.argv[3]
    if len(sys.argv) >= 5:
        frac = sys.argv[4]
    else:
        print("using default dirname """)
        frac = ""
    head = "/home/kazu/desktop/210108/Tatsumi/srlz/"+dirname+"/"
    ofile = head + "run"+runno+"united_"+frac+"spectra.pkl"
    ifile = head + "run"+runno+"united_"+frac+"monispectra.pkl"
    rfile = "./qens_hist_o_divided_by_i_"+runno+frac+denew+".pkl"
    proj = qkrodic.odata_divided_by_idata(ofile, ifile, iskde=False, denew=float(denew))
    proj.get_data()
    #proj.plot_data()
    proj.save_data(rfile)


run()