#!/usr/bin/env python
# This script balloon-estimates the re-sampled qens 1D intensity profiles with
# the kernel band widths optimized on the original data, and fits the estimated
# target density with the estimated device density.
# Kazuyoshi TATSUMI 2023/02/23
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
from qens_balloon_resample_classm import qens_balloon_resamples as qbr
#from qens_balloon_resample_org_classm_class import qens_org_classm as qoc
np.set_printoptions(suppress=True)


class qens_org_classm(qbr):
    def __init__(self, qidx, qsize, elim=[-0.03, 0.07], rsmodifier="org",
                 ishist=True, orgmodifier="org", variables=[0.655, 0.0129,
                 0.200, 0.00208], prefixes=["./", "./"], num=6400, M=160,
                 winparam=1, WinFunc='Boxcar', quiet=True, io=False,
                 runNos=[6206, 6204], ispltchk=False, isnovariablebw=False):
        self.qsize = qsize
        Nb = 1
        if ishist:
            self.outfile = "outhist" + str(runNos[0]) + "m.pkl"
        else:
            if io:
                if isnovariablebw:
                    self.outfile = "outkdeionovbw" + str(runNos[0]) + "m.pkl"
                else:
                    self.outfile = "outkdeio" + str(runNos[0]) + "m.pkl"
            else:
                if isnovariablebw:
                    if prefixes[0] == prefixes[1]:
                        self.outfile = "outkdenovbw" + str(runNos[0]) + "m.pkl"
                    else:
                        self.outfile = "outkdenovbwprefixes" + str(runNos[0]) + "m.pkl"
                else:
                    self.outfile = "outkde" + str(runNos[0]) + "m.pkl"
        qbr.__init__(self, qidx, runNos=runNos, elim=elim, Nb=Nb,
                     ishist=ishist, num=num, M=M, WinFunc=WinFunc,
                     winparam=winparam, rsmodifier=rsmodifier,
                     orgmodifier=orgmodifier, prefixes=prefixes,
                     variables=variables, quiet=quiet,
                     ispltchk=ispltchk, isnovariablebw=isnovariablebw)

    def run_for_mqs(self):
        for qidx in range(self.qsize):
            print("CHECK qidx:", qidx)
            self.qidx = qidx + 9
            if "kdeio" in self.outfile:
                self.run_eachkde_io()
            else:
                self.run_eachkde()
            if qidx == 0:
                outall = self.outall
            else:
                outall = np.concatenate((outall, self.outall), axis=0)
        if self.rank == 0:
            self.outall = outall
            self.savefile()

    def run_eachkde(self):
        self.outall = []
        for inb in range(self.Nb):
            self.kys = [self.CalcBandW(orgfile, inb=inb) for orgfile in
                        self.orgfiles]
            if self.rank == 0 and self.qidx >= 12:
                import matplotlib.pyplot as plt
                plt.plot(self.kys[0][1], self.kys[0][0])
                plt.show()
            self.DoQfandKDE(inb)
        self.outall = np.array(self.outall)
        if self.Nb > 1 and self.rank == 0:
            self.output()

    def DoQfandKDE(self, inb):
        print("ENTERING DoAfandKDE")
        xt, yt, yth = self.eachrunno(0, inb)
        xd, yd, ydh = self.eachrunno(1, inb)
        self.icorr()
        xtl, ytl = self.limit2(xt, yt, self.elim)
        xdl, ydl = self.limit2(xd, yd, self.elim)
        if inb == 0 and self.rank == 0:
            if np.sum(xtl - xdl) > 0.000001:
                print('WARNING, check x_tf - x_df')
        ydlc, ytlc = self.correction(xtl, ydl, ytl)
        #self.bg = 0.
        self.check_out(inb, self.optimize(xdl, ydlc, ytlc,
                                          variables=self.variables))
        [alpha, gamma, delta, base] = self.outall[-1][0:4]
        yqens = alpha*self.convloreorg(ydlc, gamma, xdl)
        y = yqens + delta*ydl + base
        #_yqens, _y = self.decorrection(xtl, yqens, y)
        #_ydlc, _ytlc = self.decorrection(xtl, ydlc, ytlc)
        self.ml = y
        self.yd = ydlc
        self.alpha = 0.03375
        np.random.seed(100)
        ml_data = np.zeros_like(self.ml)
        if self.rank == 0:
            yd_data, ml_data = self.generate_data()
        #yd_data, ml_data = self.correction(xtl, yd_data, ml_data)
        from mpi4py import MPI
        ml_data = MPI.COMM_WORLD.bcast(ml_data)
        self.kde(xtl, ml_data, M=self.M, winparam=self.winparam, num=self.num,
                 WinFunc=self.WinFunc, isnovariablebw=False)
        if self.rank == 0:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(2, 1)
            ax[0].plot(xtl, self.ml/np.max(self.ml))
            ax[0].plot(self.y[1], self.y[0]/np.max(self.y[0]))
            ax[1].plot(xtl, ml_data, marker='.', linewidth=0, ms=0.4)
            plt.savefig('kde_check.png')
            plt.show()



def run():
    if len(sys.argv) >= 2:
        runNos = [int(sys.argv[1]), 3540]
    else:
        runNos = [4174, 3540]
    ishist = False
    qidx = 0
    elim = [-0.025, 0.100]
    qsize = 12
    qsize = 1
    num = 8000
    M = 80
    winparam = 5
    WinFunc = 'Boxcar'
    ispltchk = True
    isnovariablebw = False
    prj = qens_org_classm(qidx, qsize, elim=elim, runNos=runNos, ishist=ishist,
                          num=num, M=M, winparam=winparam, WinFunc=WinFunc,
                          ispltchk=ispltchk, isnovariablebw=isnovariablebw)
    prj.run_for_mqs()


run()
