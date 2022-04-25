#!/usr/bin/env python
import os
import numpy as np
import pickle
import sys
sys.path.append("/home/kazu/ktpro")
from qens_fit_class_kde_idata import runkdeidata as rki


class runhistwithidata(rki):
    def __init__(self, devf, tf, idevf, itf, outfile, irwdf, alpha, elim,
                 elimw, numcycle=100, leastsq=True):
        super().__init__(devf, tf, idevf, itf, outfile, alpha, elim, elimw,
                         numcycle=numcycle, leastsq=leastsq)
        self.get_ixrw(irwdf)

    def cycle(self):
        self.outall = []
        for cyidx in range(0, self.numcycle):
            simd, simt = self.generate_data()
            tin_real, simdr = self.rebin_generated_samples(self.x, simd,
                                                           num=480)
            tin_real, simtr = self.rebin_generated_samples(self.x, simt,
                                                           num=480)
            if cyidx == 0:
                self.iyd = np.interp(self.ixrw, self.ixd, self.iyd)
                self.iyt = np.interp(self.ixrw, self.ixt, self.iyt)
                dx = tin_real[1] - tin_real[0]
            simid, simit = self.generate_data(idata=True)
            simdi = simdr/np.interp(tin_real, self.ixrw, simid)
            simti = simtr/np.interp(tin_real, self.ixrw, simit)
            simdc, simtc = self.correction(tin_real, simdi, simti)
            simdc = simdc / np.sum(simdc) / dx
            simtc = simtc / np.sum(simtc) / dx * 100.
            _out = self.optimize(tin_real, simdc, simtc,
                                 variables=[4.8e+01, 2.65e-02,
                                            3.2e+01, 7.0e-03,
                                            1.9e+01, 1.1e+01])
            self.check_out(cyidx, _out)

#    def multii(self, itf, idf, ml, yd):
#        xit, yit = self.get_idata(itf)
#        xitl, yitl = self.limit(xit, yit, self.elimw)
#        yit_ip = np.interp(self.x, xitl, yitl)
#        xid, yid = self.get_idata(idf)
#        xidl, yidl = self.limit(xid, yid, self.elimw)
#        yid_ip = np.interp(self.x, xidl, yidl)
#        return ml*yit_ip, yd*yid_ip
#
#    def get_yirw(self, itf, idf):
#        xit, yit = self.get_irwdata(itf)
#        xitl, yitl = self.limit(xit, yit, self.elimw)
#        xid, yid = self.get_irwdata(idf)
#        xidl, yidl = self.limit(xid, yid, self.elimw)
#        return np.interp(self.x, xitl, yitl), np.interp(self.x, xidl, yidl)
#
    def get_ixrw(self, ifile):
        ixrw, iyrw = self.get_irwdata(ifile)
        mergin = 0.001
        xlim = [-0.10025 - mergin, 0.14975 + mergin]
        # self.ixrw, iyrw = self.limit(ixrw, iyrw, self.elimw)
        self.ixrw, iyrw = self.limit(ixrw, iyrw, xlim)

    def get_irwdata(self, ifile):
        dataset = self.read_pkl(ifile)
        return dataset['spectra'][0, 0, :] - 2.085, dataset['spectra'][0, 1, :]

    def read_pkl(self, pklfile):
        with open(pklfile, 'rb') as f:
            dataset = pickle.load(f, encoding='latin1')
        return(dataset)


def testrun():
    np.set_printoptions(linewidth=120)
    outfile = "./outhistwidata.pkl"
    alpha = 0.5 
    devf = "./qens_kde_o_divided_by_i_6204.pkl"
    tf = "./qens_kde_o_divided_by_i_6202.pkl"
    prefix = "/home/kazu/desktop/210108/Tatsumi/pickles/0000001io/"
    idevf = prefix + "qens_kde_results_on_idata_6204.pkl"
    itf = prefix + "qens_kde_results_on_idata_6202.pkl"
    idfrw = prefix + "../../srlz/0000001io/run6204united_monispectra.pkl"
    elim = [-0.03, 0.07]
    elimw = [-0.04, 0.08]
    numcycle = 300
    binwidth1 = 0.0016
    binwidth2 = 0.0016
    if os.path.isfile(outfile):
        proj = runhistwithidata(devf, tf, idevf, itf, outfile, idfrw, alpha,
                                elim, elimw, leastsq=True, numcycle=numcycle)
        proj.loadfile()
        proj.output()
        proj.plot_distribution(binwidth1, binwidth2)
    else:
        np.random.seed(314)
        proj = runhistwithidata(devf, tf, idevf, itf, outfile, idfrw, alpha,
                                elim, elimw, leastsq=True, numcycle=numcycle)
        proj.get_xmlyd()
        proj.cycle()
        proj.output()
        proj.savefile()



#testrun()