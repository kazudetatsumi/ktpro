#!/usr/bin/env python
# This script balloon-estimates the re-sampled qens 1D intensity profiles with
# the kernel band widths optimized on the original data, and fits the estimated
# target density with the estimated device density.
# Kazuyoshi TATSUMI 2023/02/23
#import os
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
#from qens_kde_results_odata_divided_by_idata_class\
#    import odata_divided_by_idata as odbi
#from get_qlist_nova_class import get_qlist as gq
from get_resampled_data_class import Sget_qlist as sgq
from qens_kde_resampled import qens_kde_resampled as qkr
#pwd = os.getcwd()
from qens_fit_class import qens_fit as qf
#import matplotlib
#matplotlib.use('Agg')


#class qens_balloon_resamples(sgq, qkr, qf):
class qens_balloon_resamples(qf):
    def __init__(self, runNos=[6202, 6204], elim=[-0.03, 0.07], Nb=1):
        self.runNos = runNos
        self.Nb = Nb
        self.gammas = np.zeros((Nb, 2))
        self.DefineFiles()
        self.elim = [-0.03, 0.07]
        #np.set_printoptions(threshold=12, linewidth=150, suppress=True)

    def getrsspectra(self, rsfile, inb=0):
        prj = sgq(pklfile=rsfile)
        prj.load_pkl(ispython2=True)
        return prj.spectrab[inb, 0, :], prj.spectrab[inb, 1, :]

    def CalcBandW(self, rsfile, inb=0):
        proj = qkr(pklfile=rsfile)
        proj.kde(proj.spectrab[inb, 0, :], proj.spectrab[inb, 2, :])
        return proj.y

    def balloon(self, ky, sy):
        bw = np.interp(sy[0], ky[1], ky[2])
        idx = sy[1].nonzero()
        sy_nz = sy[1][idx]
        t_nz = sy[0][idx]
        yv = np.zeros_like(sy[0])
        dt = min(np.diff(sy[0]))
        for k in range(yv.shape[0]):
            yv[k] = np.sum(sy_nz * dt * self.Gauss(sy[0][k]-t_nz, bw[k]))
        yv = yv * np.max(ky[0]) / np.max(yv)
        return yv

    def Gauss(self, x, w):
        y = 1 / (2 * np.pi)**2 / w * np.exp(-x**2 / 2 / w**2)
        return y

    def eachrunno(self, fidx, inb):
        sy = self.getrsspectra(self.rsfiles[fidx], inb)
        syb = self.balloon(self.kys[fidx], sy)
        return sy[0], syb, sy[1]

    def DefineFiles(self):
        self.rsfiles = []
        self.orgfiles = []
        for runno in self.runNos:
            self.rsfiles.append("./run" + str(runno) + "spectrab.pkl")
            self.orgfiles.append("./run" + str(runno) + "spectraorg.pkl")

    def DoQf(self, inb):
        qf.__init__(self, 'dummy', 'dummy', self.elim, showplot=False, leastsq=False, quiet=True)
        xt, yt, yth = self.eachrunno(0, inb)
        xd, yd, ydh = self.eachrunno(1, inb)
        self.icorr()
        self.x_tf, self.y_tf = self.limit2(xt, yt, self.elim)
        self.x_df, self.y_df = self.limit2(xd, yd, self.elim)
        self.correction()
        self.bg = 0.
        self.optimize(variables=[0.655, 0.0129, 0.200, 0.00208])

    def run(self):
        self.kys = [self.CalcBandW(orgfile, inb=0) for orgfile in self.orgfiles]
        for inb in range(self.Nb):
            self.DoQf(inb)
            self.gammas[inb, 0] = self.out[1]
            self.gammas[inb, 1] = self.out[3]
        print(np.mean(self.gammas, axis=0))
        print(np.std(self.gammas, axis=0))

 
def testrun():
    Nb = 4
    elim = [-0.03, 0.07]
    prj = qens_balloon_resamples(runNos=[6202, 6204], elim=elim, Nb=Nb)
    prj.run()


testrun()
