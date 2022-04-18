#!/usr/bin/env python
# example command: qens_qselected_using_class.py  6204 000025io 80 1 Boxcar 0875
import sys
sys.path.append("/home/kazu/ktpro")
import qens_class_fort_mpi as qc


def run():
    if len(sys.argv) >= 2:
        runno = sys.argv[1]
    datadir="./"
    WinFunc="Boxcar"
    save_file = "qens_sim_"+runno+".pkl"
    output_file = "./qens_sim_kde_"+runno+".pkl"
    figname = "qens_sim_out_"+runno+".png"
    M=160
    winparam=1
    proj = qc.qens(datadir, save_file, qsel=True, winparam=winparam, M=M,
                   WinFunc=WinFunc, figname=figname, showplot=False)
    proj.select_spectra()
    proj.add_shift_de()
    proj.run_ssvkernel()
    proj.plotter()
    proj.save_output(output_file)


run()
