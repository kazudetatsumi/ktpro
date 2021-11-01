#!/usr/bin/env python
# qens_fit_results for all used  fracs of event data.
# Kazuyoshi TATSUMI 2021.
# example: qens_fit_kde_hist_using_class.py 000025new
import sys
sys.path.append("/home/kazu/ktpro")
import qens_fit_kde_hist_class as qfkhc
argvs = sys.argv


def run():
    head = "/home/kazu/desktop/210108/Tatsumi/"
    dirname = argvs[1]
    fsts = ["0.375", "0.5", "0.625", "0.75", "0.875", "1"]
    outf = "kh_fit_results.txt"
    prj = qfkhc.sqrun_kde_hist()
    prj.fitalldata(head, dirname, fsts, outf)
    prj.plotter()


run()
