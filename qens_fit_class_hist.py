#!/usr/bin/env python
import pickle
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
from qens_fit_class_hist_noidata import runhistnoidata as rhn


class runhist(rhn):
    def __init__(self, devf, tf, idevf, itf, idfrw, itfrw, elim, elimw,
                 numcycle=100):
        self.elim = elim
        self.elimw = elimw
        self.devf = devf
        self.tf = tf
        self.idevf = idevf
        self.itf = itf
        self.idfrw = idfrw
        self.itfrw = itfrw
        self.elim = elim
        self.numcycle = numcycle

    def get_xmlyd(self):
        x, yd, yt = self.preprocess()
        out = self.optimize(x, yd, yt,
                            variables=[2.18704786e-04, 1.67980295e-02,
                                       4.92405238e-05, 1.88866588e-03,
                                       1.21127501e-01, 5.02759930e-02])
        ml = self.reconstruct(x, yd, out)
        mld, ydd = self.decorrection(x, ml, yd)
        self.x = x
        self.ml, self.yd = self.multii(self.itf, self.idevf, mld, ydd)

    def cycle(self):
        self.outall = np.zeros((self.numcycle, 6))
        yitrw, yidrw = self.get_yirw(self.itfrw, self.idfrw)
        for cyidx in range(0, self.numcycle):
            simd, simt = self.generate_data()
            out = self.optimize(self.x, simd/yidrw, simt/yitrw,
                                variables=[6.11704786e-06, 2.51980295e-02,
                                           1.55405238e-06, 4.28866588e-03,
                                           7.97127501e-03, 3.52759930e-01])
            self.outall[cyidx, :] = out
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
        print(self.outnonneg.shape[0], "/", self.numcycle)

    def decorrection(self, x, ml, yd):
        x = x + 2.085
        mld = ml * (self.k[0] + self.k[1]*x + self.k[2]*x**2 + self.k[3]*x**3)
        ydd = yd * (self.k[0] + self.k[1]*x + self.k[2]*x**2 + self.k[3]*x**3)
        return mld, ydd

    def multii(self, itf, idf, ml, yd):
        xit, yit = self.get_idata(itf)
        xitl, yitl = self.limit(xit, yit, self.elimw)
        yit_ip = np.interp(self.x, xitl, yitl)
        xid, yid = self.get_idata(idf)
        xidl, yidl = self.limit(xid, yid, self.elimw)
        yid_ip = np.interp(self.x, xidl, yidl)
        return ml*yit_ip, yd*yid_ip

    def get_yirw(self, itf, idf):
        xit, yit = self.get_irwdata(itf)
        xitl, yitl = self.limit(xit, yit, self.elimw)
        xid, yid = self.get_irwdata(idf)
        xidl, yidl = self.limit(xid, yid, self.elimw)
        return np.interp(self.x, xitl, yitl), np.interp(self.x, xidl, yidl)

    def get_irwdata(self, ifile):
        dataset = self.read_pkl(ifile)
        return dataset['spectra'][0, 0, :] - 2.085, dataset['spectra'][0, 1, :]

    def read_pkl(self, pklfile):
        with open(pklfile, 'rb') as f:
            dataset = pickle.load(f, encoding='latin1')
        return(dataset)


def testrun():
    prefix = "/home/kazu/desktop/210108/Tatsumi/"
    np.set_printoptions(linewidth=120)
    devf = "./qens_kde_o_divided_by_i_6204.pkl"
    tf = "./qens_kde_o_divided_by_i_6202.pkl"
    idevf = prefix + "pickles/0000001io/qens_kde_results_on_idata_6204.pkl"
    itf = prefix + "pickles/0000001io/qens_kde_results_on_idata_6202.pkl"
    idfrw = prefix + "srlz/0000001io/run6204united_monispectra.pkl"
    itfrw = prefix + "srlz/0000001io/run6202united_monispectra.pkl"
    elim = [-0.03, 0.07]
    elimw = [-0.06, 0.10]
    proj = runhist(devf, tf, idevf, itf, idfrw, itfrw, elim, elimw,
                   numcycle=138)
    proj.get_xmlyd()
    proj.cycle()


testrun()
