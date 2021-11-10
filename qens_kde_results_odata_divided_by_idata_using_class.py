#!/usr/bin/env python
import sys
sys.path.append("/home/kazu/ktpro")
import qens_kde_results_odata_divided_by_idata_class as qkrodic


def run():
    head = "/home/kazu/desktop/210108/Tatsumi/pickles/000025io/"
    ofile = head + "qens_run6202united_kde_results_on_data_qsel.pkl"
    ifile = head + "qens_kde_results_on_idata_6202.pkl"
    rfile = head + "qens_kde_o_divided_by_i_6202.pkl"
    proj = qkrodic.odata_divided_by_idata(ofile, ifile)
    #head = "/home/kazu/desktop/210108/Tatsumi/srlz/000025io/"
    #ofile = head + "run6202united_spectra.pkl"
    #ifile = head + "run6202_monispectra.pkl"
    #rfile = head + "qens_hist_o_divided_by_i_6202.pkl"
    #proj = qkrodic.odata_divided_by_idata(ofile, ifile, iskde=False)
    proj.get_data()
    proj.plot_data()
    proj.save_data(rfile)


run()

