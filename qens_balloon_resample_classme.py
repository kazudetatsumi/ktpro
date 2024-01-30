
# This script balloon-estimates the re-sampled qens 1D intensity profiles with
# the kernel band widths optimized on the original data, and fits the estimated
# target density with the estimated device density.
# Kazuyoshi TATSUMI 2023/02/23
from mpi4py import MPI
import numpy as np
import pickle
import sys
sys.path.append("/home/kazu/ktpro")
from qens_balloon_resample_class import Sqens_balloon_resamples as sqkr


class qens_balloon_resamples(sqkr):
    def __init__(self, qidx, runNos=[6202, 6204], elim=[-0.03, 0.07], Nb=1,
                 ishist=False, num=6400, rsmodifier="b", orgmodifier="orge",
                 prefixes=["./", "./"],
                 variables=[0.655, 0.0129, 0.200, 0.00208], quiet=False,
                 ispltchk=False):
        self.qidx = qidx
        self.runNos = runNos
        self.Nb = Nb
        self.gammas = np.zeros((Nb, 2))
        self.elim = elim
        self.ishist = ishist
        self.num = num
        self.rsmodifier = rsmodifier
        self.orgmodifier = orgmodifier
        self.prefixes = prefixes
        self.variables = variables
        self.quiet = quiet
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
        self.leastsq = False
        self.ispltchk = ispltchk
        self.DefineFiles()

    def getrsspectra(self, rsfile, inb=0):
        super(sqkr, self).__init__(pklfile=rsfile)
        print("getrsspectra: chkm slicing spectrab at qidx")
        return self.spectrab[inb, 0, self.qidx],\
            self.spectrab[inb, 1, self.qidx],\
            self.spectrab[inb, 2, self.qidx]

    def CalcBandW(self, orgfile, inb=0):
        if self.ishist:
            if self.rank == 0 and not self.quiet:
                print("skipping KDE because ishist", self.ishist)
            self.y = "dummy"
        else:
            super(sqkr, self).__init__(pklfile=orgfile)
            print("CalcBandW: chkm slicing spectrab at qidx")
            self.spectrab = self.spectrab[:, :, self.qidx, :]
            self.kde(self.spectrab[inb, 0, :], self.spectrab[inb, 2, :],
                     num=self.num)
        return self.y

    def geterrorbars(self):
        if self.ishist:
            with open("outallhst.pkl."+str(self.qidx), 'rb') as f:
                dat = pickle.load(f)['out']
        else:
            with open("outallkde.pkl."+str(self.qidx), 'rb') as f:
                dat = pickle.load(f)['out']
        return np.std(dat, axis=0)

    def geterrorbars_io(self):
        with open("outallkdeio.pkl."+str(self.qidx), 'rb') as f:
            data = pickle.load(f)
            if 'x' in data.keys():
                self.xe = data['x']
        return np.std(data['out'], axis=0)

    def res(self, coeffs, x, d, t):
        if len(coeffs) == 6:
            [alpha1, gamma1, alpha2, gamma2,  delta, base] = coeffs
            y = alpha1*self.convlore(d, gamma1, x)\
                + alpha2*self.convlore(d, gamma2, x)\
                + delta*d + base
        if len(coeffs) == 5:
            [alpha1, gamma1, alpha2, gamma2,  delta] = coeffs
            y = alpha1*self.convlore(d, gamma1, x)\
                + alpha2*self.convlore(d, gamma2, x)\
                + delta*d + self.bg
        if len(coeffs) == 4:
            [alpha, gamma, delta, base] = coeffs
            y = alpha*self.convlore(d, gamma, x)\
                + delta*d + base
        if len(coeffs) == 3:
            [alpha, gamma, delta] = coeffs
            y = alpha*self.convlore(d, gamma, x)\
                + delta*d + self.bg
        #xl, dif = self.limit2(x, (t-y)/self.etl, self.elim)
        xl, dif = self.limit2(x, t-y, self.elim)
        xl, e = self.limit2(x, self.etl, self.elim)
        dif = dif[np.nonzero(e)]/e[np.nonzero(e)]
        #xl, dif = self.limit2(x, t-y, self.elim)
        return dif

    def run_eachkde(self):
        self.etl = self.geterrorbars()
        self.outall = []
        for inb in range(self.Nb):
            self.kys = [self.CalcBandW(orgfile, inb=inb) for orgfile in
                        self.orgfiles]
            self.DoQf(inb)
        self.outall = np.array(self.outall)
        if self.Nb > 1 and self.rank == 0:
            self.output()

    def run_eachkde_io(self):
        self.etl = self.geterrorbars_io()
        self.outall = []
        self.kyos = [self.CalcBandW(orgfile, inb=0) for orgfile in
                     self.orgfiles]
        self.check_idata()
        self.kyios = [self.io(kyo, kyi) for kyo, kyi in
                      zip(self.kyos, self.kyis)]
        print("CHK2", self.kyios[0].shape, self.kyios[1].shape)
        self.DoQfio(0)
        self.outall = np.array(self.outall)
        if self.Nb > 1 and self.rank == 0:
            self.output()


def testrun():
    np.set_printoptions(suppress=True)
    Nb = 1
    elim = [-0.03, 0.07]
    ishist = True
    num = 6400
    rsmodifier = "org"
    variables = [0.8, 0.01, 0.24, 0.0002, 0.001, 1.2]
    variables = [0.655, 0.0129, 0.200, 0.00208]
    preprefix = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/"
    #prefix = preprefix + "100/"
    #prefix = preprefix + "0125/back/test5/6208/nonmpi_test/"
    #prefix = preprefix + "0125/back/test5/6208/whole/"
    prefix = preprefix + "test/0125/back/test5/6208/nonmpi_test/dnapca03/"
    #prj = qens_balloon_resamples(runNos=[6202, 6204], elim=elim, Nb=Nb,
    prj = qens_balloon_resamples(runNos=[6207, 6204], elim=elim, Nb=Nb,
                                 ishist=ishist, num=num, rsmodifier=rsmodifier,
                                 prefixes=[prefix, prefix], variables=variables)
    #print(qens_balloon_resamples.__mro__)
    prj.run()
    prj.ishist = False
    prj.run()


#testrun()
