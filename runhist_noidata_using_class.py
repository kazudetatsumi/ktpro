#!/usr/bin/env python
# This script two-lorentzian-fits kde results, and reconstructs the target
# function with the fitting parameter.
# The reconstructed ML function will be used to generate simulated data.
# example: qens_fit_kde_vs_hist_2lore_reconst_using_class.py 6202 6204 0875
# Kazuyoshi TATSUMI 2022/04/16
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
from qens_fit_class import qens_fit as qf
#from qens_calss_fort_mpi import qens as qens


class runhist(qf):
    def __init__(self, devf, tf, elim, elimr, numcycle=100):
        self.elimo = elim
        self.elimr = elimr
        self.devfo = devf
        self.tfo = tf
        self.numcycle = numcycle
        qf.__init__(self, devf, tf, self.elimo, showplot=False)
        self.icorr()
        self.preprocess(doicorr=True)
        self.optimize(variables=[2.18704786e-04, 1.67980295e-02,
                                 4.92405238e-05, 1.88866588e-03,
                                 1.21127501e-01, 5.02759930e-02],
                      figname="qens_kde_fit2.png")
        self.outo = self.out

    def cycle(self):
        self.outall = np.zeros((self.numcycle, 6))
        for cyidx in range(0, self.numcycle):
            self.devf = self.devfo
            self.tf = self.tfo
            self.elim = self.elimr
            self.outall[cyidx, :] = self.out
            self.out = self.outo
            self.reconstruct(elim=self.elimr, check=False)
            self.generate_data(None, None, check=False)
            self.devf = "./qens_sim_6204.pkl"
            self.tf = "./qens_sim_6202.pkl"
            self.elim = self.elimo
            self.preprocessh(doicorr=False)
            self.optimize(variables=[6.11704786e-06, 2.51980295e-02,
                                     1.55405238e-06, 4.28866588e-03,
                                     7.97127501e-03, 3.52759930e-01],
                          figname="qens_kde_fit2.png")
        print(np.average(self.outall[:, 1]), np.std(self.outall[:, 1]))
        print(np.average(self.outall[:, 3]), np.std(self.outall[:, 3]))
        mask = np.where((self.outall[:, 0] > 0) & (self.outall[:, 1] > 0)
                        & (self.outall[:, 2] > 0) & (self.outall[:, 3] > 0)
                        & (self.outall[:, 4] > 0) & (self.outall[:, 5] > 0))
        self.outnonneg = self.outall[mask]
        print(np.average(self.outnonneg[:, 1]), "+/-",
              np.std(self.outnonneg[:, 1]))
        print(np.average(self.outnonneg[:, 3]), "+/-",
              np.std(self.outnonneg[:, 3]))


def testrun():
    np.set_printoptions(linewidth=120)
    devf = "./qens_kde_o_divided_by_i_6204.pkl"
    tf = "./qens_kde_o_divided_by_i_6202.pkl"
    elim = [-0.03, 0.07]
    elimr = [-0.06, 0.10]
    proj = runhist(devf, tf, elim, elimr, numcycle=640)
    proj.cycle()


testrun()
