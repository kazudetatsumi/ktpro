#!/usr/bin/env python
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from mpi4py import MPI
sys.path.append("/home/kazu/ktpro")
#from qens_class_fort_mpi import qens as qc
from qens_fit_class_hist_widata import runhistwithidata as rhw


#class runkdewithidata(rhw, qc):
class runkdewithidata(rhw):
    def __init__(self, devf, tf, idevf, itf, outfile, irwdf, alpha, elim,
                 elimw, numcycle=100, leastsq=True):
        super().__init__(devf, tf, idevf, itf, outfile, irwdf, alpha, elim,
                         elimw, numcycle=numcycle, leastsq=leastsq)

    def cycle(self):
        self.outall = []
        for cyidx in range(0, self.numcycle):
            # odata
            simt = np.zeros(self.ml.shape)
            simd = np.zeros(self.yd.shape)
            if self.rank == 0:
                simd, simt = self.generate_data()
            simt = MPI.COMM_WORLD.bcast(simt)
            simd = MPI.COMM_WORLD.bcast(simd)
            # idata
            if cyidx == 0:
                self.iyd = np.interp(self.ixrw, self.ixd, self.iyd)
                self.iyt = np.interp(self.ixrw, self.ixt, self.iyt)
            simit = np.zeros(self.iyt.shape)
            simid = np.zeros(self.iyd.shape)
            if self.rank == 0:
                simid, simit = self.generate_data(idata=True)
            simit = MPI.COMM_WORLD.bcast(simit)
            simid = MPI.COMM_WORLD.bcast(simid)
            # kde on idata
            self.kde(self.ixrw, simid)
            simiyd = self.y[0]
            self.kde(self.ixrw, simit)
            simiyt = self.y[0]
            ix = self.y[1]
            # kde on odata
            self.kde(self.x, simd)
            self.dt = self.y[1][1]-self.y[1][0]
            simyd = self.y[0]
            simydi = simyd/np.interp(self.y[1], ix, simiyd)
            self.kde(self.x, simt)
            simyt = self.y[0]
            simyti = simyt/np.interp(self.y[1], ix, simiyt)
            simydc, simytc = self.correction(self.y[1], simydi, simyti)
            simydc = simydc/np.sum(simydc)/self.dt
            simytc = simytc/np.sum(simytc)/self.dt*100.
            _out = self.optimize(self.y[1], simydc, simytc,
                                 variables=[4.8e+01, 2.65e-02,
                                            3.2e+01, 7.0e-03,
                                            1.9e+01, 1.1e+01])
            self.check_out(cyidx, _out)


def testrun():
    np.set_printoptions(linewidth=120)
    outfile = "./outkdewidata.pkl"
    alpha = 0.5
    devf = "./qens_kde_o_divided_by_i_6204.pkl"
    tf = "./qens_kde_o_divided_by_i_6202.pkl"
    prefix = "/home/kazu/desktop/210108/Tatsumi/pickles/0000001io/"
    idevf = prefix + "qens_kde_results_on_idata_6204.pkl"
    itf = prefix + "qens_kde_results_on_idata_6202.pkl"
    idfrw = prefix + "../../srlz/0000001io/run6204united_monispectra.pkl"
    elim = [-0.03, 0.07]
    elimw = [-0.04, 0.08]
    numcycle = 3
    binwidth1 = 0.0016
    binwidth2 = 0.0016
    if os.path.isfile(outfile):
        proj = runkdewithidata(devf, tf, idevf, itf, outfile, idfrw, alpha,
                               elim, elimw, leastsq=True, numcycle=numcycle)
        if proj.rank == 0:
            proj.loadfile()
            proj.output()
            proj.plot_distribution(binwidth1, binwidth2)
    else:
        np.random.seed(314)
        proj = runkdewithidata(devf, tf, idevf, itf, outfile, idfrw, alpha,
                               elim, elimw, leastsq=True, numcycle=numcycle)
        proj.get_xmlyd()
        proj.cycle()
        if proj.rank == 0:
            proj.output()
            proj.savefile()


#testrun()
