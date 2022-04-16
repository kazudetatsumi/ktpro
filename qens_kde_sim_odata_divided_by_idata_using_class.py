#!/usr/bin/env python
import sys
sys.path.append("/home/kazu/ktpro")
import qens_kde_results_odata_divided_by_idata_class as qkrodic

# example command:  qens_kde_results_odata_divided_by_idata_using_class.py 6204 000025io 0875

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
        frac = sys.argv[3]
    else:
        print("using default frac """)
        frac = ""
    headorg = "/home/kazu/desktop/210108/Tatsumi/pickles/"+dirname+"/"
    head = "./"
    ofile = head + "qens_sim_kde_"+runno+frac+".pkl"
    ifile = headorg + "qens_kde_results_on_idata_"+runno+frac+".pkl"
    rfile = head + "qens_kde_sim_o_divided_by_i_"+runno+frac+".pkl"
    proj = qkrodic.odata_divided_by_idata(ofile, ifile)
    proj.get_data()
    #proj.plot_data()
    proj.save_data(rfile)


run()
