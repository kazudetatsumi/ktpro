#!/usr/bin/env python
# This script balloon-estimates the re-sampled qens 1D intensity profiles with
# the kernel band widths optimized on the original data, and fits the estimated
# target density with the estimated device density.
# Kazuyoshi TATSUMI 2023/02/23
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
from mpi4py import MPI
from qens_kde_resampled import qens_kde_resampled as qkr


class Sqens_balloon_resamples(qkr):
    def __init__(self, runNos=[6202, 6204], elim=[-0.03, 0.07], Nb=1,
                 ishist=False, num=6400, rsmodifier="b", orgmodifier="org",
                 prefixes=["./", "./"], variables=[0.655, 0.0129, 0.200,
                 0.00208], quiet=False, ispltchk=False, isnovariablebw=False):
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
        self.isnovariablebw = isnovariablebw
        self.DefineFiles()

    def getrsspectra(self, rsfile, inb=0):
        qkr.__init__(self, pklfile=rsfile)
        return self.spectrab[inb, 0, :], self.spectrab[inb, 1, :],\
            self.spectrab[inb, 2, :]

    def CalcBandW(self, orgfile, inb=0):
        if self.ishist:
            if self.rank == 0 and not self.quiet:
                print("skipping KDE because ishist", self.ishist)
            self.y = "dummy"
        else:
            qkr.__init__(self, pklfile=orgfile)
            print(orgfile)
            self.kde(self.spectrab[inb, 0, :], self.spectrab[inb, 2, :],
                     num=self.num, M=self.M, winparam=self.winparam,
                     isnovariablebw=self.isnovariablebw)
        return self.y

    def balloon(self, ky, sy):
        # if ky[2] is np.array, then ky is the result of vKDE.
        # else, ky is the result of the bandwidth of the KDE is a constant over
        # the whole energies. In the latter case, the bandwidth ndarray is
        # preapared to be suit for the balloon estimation.
        if isinstance(ky[2], np.ndarray):
            print("balloon estimation for variable bw")
            bw = np.interp(sy[0], ky[1], ky[2])
        else:
            print("balloon estimation for non-variable bw")
            bw = np.ones(sy[0].shape[0])*ky[2]
        idx = sy[1].nonzero()
        sy_nz = sy[1][idx]
        t_nz = sy[0][idx]
        yv = np.zeros_like(sy[0])
        dt = min(np.diff(sy[0]))
        for k in range(yv.shape[0]):
            yv[k] = np.sum(sy_nz * dt * self.Gauss(sy[0][k]-t_nz, bw[k]))
        ##yv = yv * np.max(ky[0]) / np.max(yv)
        return yv

    def Gauss(self, x, w):
        return 1.0/((2.0*np.pi)**0.5*w)*np.exp(-(x/w)**2/2.0)

    def DefineFiles(self):
        self.rsfiles = []
        self.orgfiles = []
        for runno, prefix in zip(self.runNos, self.prefixes):
            self.rsfiles.append(prefix + "run" + str(runno) + "spectra" +
                                self.rsmodifier + ".pkl")
            self.orgfiles.append(prefix + "run" + str(runno) + "spectra" +
                                 self.orgmodifier + ".pkl")
        if self.rank == 0:
            print("rsfiles:", self.rsfiles)
            print("orgfiles:", self.orgfiles)

    def eachrunno(self, fidx, inb):
        sy = self.getrsspectra(self.rsfiles[fidx], inb)
        if self.ishist:
            return sy[0], sy[1], sy[1]
        else:
            syb = self.balloon(self.kys[fidx], sy)
            return sy[0], syb, sy[1]

    def eachrunno_io(self, fidx, inb):
        sy = self.getrsspectra(self.rsfiles[fidx], inb)
        syb = self.balloon(self.kys[fidx], [sy[0], sy[2], sy[1]])
        return syb, sy[0]

    def DoQf(self, inb):
        xt, yt, yth = self.eachrunno(0, inb)
        xd, yd, ydh = self.eachrunno(1, inb)
        self.icorr()
        xtl, ytl = self.limit2(xt, yt, self.elim)
        xdl, ydl = self.limit2(xd, yd, self.elim)
        if inb == 0 and self.rank == 0:
            if np.sum(xtl - xdl) > 0.000001:
                print('WARNING, check x_tf - x_df')
        ydlc, ytlc = self.correction(xtl, ydl, ytl)
        self.bg = 0.
        self.check_out(inb, self.optimize(xdl, ydlc, ytlc,
                                          variables=self.variables))
        if self.rank == 0 and self.ispltchk:
            import matplotlib.pyplot as plt
            [alpha, gamma, delta, base] = self.outall[-1][0:4]
            yqens = alpha*self.convloreorg(ydlc, gamma, xdl)
            y = yqens + delta*ydl + base
            plt.plot(xdl*1000, y, c='k')
            plt.plot(xdl*1000, ytlc, c='b', label='ytlc@qidx'+str(self.qidx))
            plt.plot(xdl*1000, yqens, ls='dotted', c='k')
            plt.ylabel('Intensity (Arb. Units)')
            plt.xlabel(r'$Energy\ (\mu eV)$')
            plt.legend()
            plt.show()
        #if self.rank == 0:
        #    import matplotlib.pyplot as plt
        #    [alpha, gamma, delta, base] = self.outall[-1][0:4]
        #    yqens = alpha*self.convloreorg(ydlc, gamma, xdl)
        #    y = yqens + delta*ydl + base
        #    plt.plot(xdl*1000, y, c='k')
        ##    #plt.scatter(xdl*1000, ytlc, marker='o', s=18, fc='None', ec='k')
        #    plt.plot(xdl*1000, ytlc, c='b')
        #    plt.plot(xdl*1000, yqens, ls='dotted', c='k')
        #    plt.ylabel('Intensity (Arb. Units)')
        #    plt.xlabel(r'$Energy\ (\mu eV)$')
        #    plt.show()

    def DoQfio(self, inb):
        self.icorr()
        print("CHKKK", self.elim)
        xtl, ytl = self.limit2(self.kyos[0][1], self.kyios[0], self.elim)
        xdl, ydl = self.limit2(self.kyos[1][1], self.kyios[1], self.elim)
        print('CHK', ytl.shape)
        if inb == 0 and self.rank == 0:
            if np.sum(xtl - xdl) > 0.000001:
                print('WARNING, check x_tf - x_df')
        ydlc, ytlc = self.correction(xtl, ydl, ytl)
        de = xtl[1] - xtl[0]
        #tmpint = np.sum(ydlc)
        #ydlc /= tmpint * de
        #ytlc /= tmpint * de
        ydlc *= 1000000.
        ytlc *= 1000000.
        self.bg = 0.
        try:
            if self.xe.shape[0] != xtl.shape[0]:
                print("interpolating errors")
                self.etl = np.interp(xtl, self.xe, self.etl)
        except AttributeError:
            print('self.xe does not exist.')
        self.check_out(inb, self.optimize(xdl, ydlc, ytlc,
                                          variables=self.variables))

    def CI_of_intensities(self):
        if not self.ishist:
            self.kys = [self.CalcBandW(self.orgfiles[0], inb=0)]
        ytlcs = []
        for inb in range(self.rank*self.Nb//self.size,
                         (self.rank+1)*self.Nb//self.size):
            xt, yt, yth = self.eachrunno(0, inb)
            self.icorr()
            xtl, ytl = self.limit2(xt, yt, self.elim)
            dummy, ytlc = self.correction(xtl, ytl, ytl)
            ytlcs.append(ytlc)
        outallt = np.array(self.comm.gather(np.array(ytlcs).flatten(),
                           root=0))
        if self.rank == 0:
            self.outall = outallt.reshape((self.Nb, -1))

    def CI_of_intensities_io(self):
        #self.kys = [self.CalcBandW(orgfile, inb=0) for orgfile in self.orgfiles
        #            ]
        self.kys = [self.CalcBandW(self.orgfiles[0], inb=0)]
        self.check_idata()
        ytlcs = []
        for inb in range(self.rank*self.Nb//self.size,
                         (self.rank+1)*self.Nb//self.size):
            self.kyos = [self.eachrunno_io(fidx, inb) for fidx in range(1)]
            self.kyios = [self.io(kyo, kyi) for kyo, kyi in
                          zip(self.kyos, self.kyis)]
            self.icorr()
            xtl, ytl = self.limit2(self.kyos[0][1], self.kyios[0], self.elim)
            dummy, ytlc = self.correction(xtl, ytl, ytl)
            ytlcs.append(ytlc)
        outallt = np.array(self.comm.gather(np.array(ytlcs).flatten(), root=0))
        self.x = xtl
        if self.rank == 0:
            self.outall = outallt.reshape((self.Nb, -1))

    def check_idata(self):
        print("chek_idata")
        self.odata = False
        self.kyis = []
        for rsfile in self.rsfiles:
            print(rsfile)
            if ".pkl." in rsfile:
                self.pklfile = rsfile.split(".pkl.")[0]+".pkl.moni"
            else:
                self.pklfile = rsfile + ".moni"
            print(self.pklfile)
            self.read_pkl()
            self.select_spectra()
            for i in range(10):
                # the parameters of kde, M, num, winparam can be altered
                # by setting the corresponding options in the argument.
                # M=5 and num=800 seem preferable for smooth energy
                # distributions of the monitored neutron beam. But, I do not
                # change the default value at the present time (2023/07/20)
                self.kde(self.selected_energy, self.selected_spectra, M=160,
                         winparam=1, num=self.num)
                # if self.rank == 0:
                #     import matplotlib.pyplot as plt
                #     fig, ax = plt.subplots()
                #     ax2 = ax.twinx()
                #     ax.plot(self.y[1], self.y[0])
                #     ax2.plot(self.y[1], self.y[2])
                #     plt.show()
            self.kyis.append(self.y)

    def io(self, kyo, kyi):
        if self.rank == 0:
            print("CHK NONZERO", np.where(np.interp(kyo[1], kyi[1], kyi[0])
                                          == 0.))
        return kyo[0]/np.interp(kyo[1], kyi[1], kyi[0])

    def run(self):
        self.kys = [self.CalcBandW(orgfile, inb=0) for orgfile in self.orgfiles
                    ]
        self.outall = []
        for inb in range(self.rank*self.Nb//self.size,
                         (self.rank+1)*self.Nb//self.size):
            self.DoQf(inb)
        ##MPI Gathering the result from each rank
        outallt = np.zeros((self.size*np.array(self.outall).size))
        self.comm.Allgather(np.array(self.outall).flatten(), outallt)
        self.outall = outallt.reshape((self.Nb, -1))
        if self.Nb > 1 and self.rank == 0:
            self.output()

    def run_io(self):
        self.kys = [self.CalcBandW(orgfile, inb=0) for orgfile in self.orgfiles
                    ]
        self.check_idata()
        self.outall = []
        for inb in range(self.rank*self.Nb//self.size,
                         (self.rank+1)*self.Nb//self.size):
            self.kyos = [self.eachrunno_io(fidx, inb) for fidx in range(2)]
            self.kyios = [self.io(kyo, kyi) for kyo, kyi in
                          zip(self.kyos, self.kyis)]
            self.DoQfio(inb)
        ##MPI Gathering the result from each rank
        outallt = np.zeros((self.size*np.array(self.outall).size))
        self.comm.Allgather(np.array(self.outall).flatten(), outallt)
        self.outall = outallt.reshape((self.Nb, -1))
        if self.Nb > 1 and self.rank == 0:
            self.output()

    def run_eachkde(self):
        self.outall = []
        for inb in range(self.Nb):
            self.kys = [self.CalcBandW(orgfile, inb=inb) for orgfile in
                        self.orgfiles]
            self.DoQf(inb)
        self.outall = np.array(self.outall)
        if self.Nb > 1 and self.rank == 0:
            self.output()

    def run_eachkde_io(self):
        self.outall = []
        self.kyos = [self.CalcBandW(orgfile, inb=0) for orgfile in
                     self.orgfiles]
        self.check_idata()
        self.kyios = [self.io(kyo, kyi) for kyo, kyi in
                      zip(self.kyos, self.kyis)]
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
                                 prefixes=[prefix, prefix],
                                 variables=variables)
    #print(qens_balloon_resamples.__mro__)
    prj.run()
    prj.ishist = False
    prj.run()


#testrun()
