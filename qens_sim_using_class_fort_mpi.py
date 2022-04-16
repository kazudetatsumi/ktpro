#!/usr/bin/env python
# example command: qens_qselected_using_class.py  6204 000025io 80 1 Boxcar 0875
import sys
sys.path.append("/home/kazu/ktpro")
import qens_class_fort_mpi as qc


def run():
    datadir="./"
    WinFunc="Boxcar"
    save_file = "qens_sim_6204.pkl"
    output_file = "./qens_sim_kde_6204.pkl"
    figname = "qens_sim_out_6204.png"
    M=160
    winparam=10
    proj = qc.qens(datadir, save_file, qsel=True, winparam=winparam, M=M,
                   WinFunc=WinFunc, figname=figname, showplot=True)
    proj.select_spectra()
    #proj.add_shift_de()
    proj.run_ssvkernel()
    proj.plotter()
    proj.save_output(output_file)


run()
