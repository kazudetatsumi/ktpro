#i!/usr/bin/env python
import numpy as np
import pickle
from mpi4py import MPI
import sys
sys.path.append("/home/kazu/ktpro")
from change_h_pos_in_poscar_class import change_hpos as ch


class change_hpos(ch):
    def __init__(self, infile, std, edgelength, nx, enefile=None,
                 shift=[0., 0., 0.], hshift=[0., 0., 0.], rg=None, prim=True,
                 oldedgelength=False, deuterium=False):
        super(change_hpos, self).__init__(infile, std, edgelength, nx,
                                          enefile=enefile, shift=shift,
                                          hshift=hshift, rg=rg, prim=prim,
                                          oldedgelength=oldedgelength,
                                          deuterium=deuterium)

    def GetEigen(self, Issave=False):
        rank = MPI.COMM_WORLD.Get_rank()
        self.E, self.U = np.linalg.eigh(self.H)
        self.E *= self.Eh * 1000.
        if rank == 0:
            print(self.E[0:15] - np.min(self.E))
            if Issave:
                dataset = {}
                dataset['E'] = self.E
                dataset['U'] = self.U
                with open("./save_eigen.pkl", 'wb') as f:
                    pickle.dump(dataset, f, 4)

    def Integrateoverallsolidangle2(self, qlength, Iswrite=True):
        comm = MPI.COMM_WORLD
        rank = MPI.COMM_WORLD.Get_rank()
        size = MPI.COMM_WORLD.Get_size()
        nmesh = 16
        coss, weights = np.polynomial.legendre.leggauss(nmesh)
        thetas = np.arccos(coss)
        phis = 2.0*np.pi*np.arange(0, nmesh*2)/nmesh*2
        for it, theta in enumerate(thetas):
            #for ip, phi in enumerate(phis):
            phi = phis[rank]
            print("theta number:", it, "my rank:", rank)
            # domega = np.sin(theta)*dtheta*dphi
            # q = qlength * np.array([np.sin(theta)*np.cos(phi),
            #                         np.sin(theta)*np.sin(phi),
            #                         np.cos(theta)])
            q = qlength * np.array([np.cos(theta), np.sin(theta)*np.cos(phi),
                                    np.sin(theta)*np.sin(phi)])
            self.GetTransitionMatrix(q, Isplot=False, Iscalledbyigos=True)
            if it == 0:
                IntegratedSqw = np.zeros((self.sqw.shape[0]))
                IntegratedSqwg = np.zeros((self.sqw.shape[0]*size))
            IntegratedSqw += self.sqw*weights[it]
        comm.Allgather([IntegratedSqw, MPI.DOUBLE],
                       [IntegratedSqwg, MPI.DOUBLE])
        IntegratedSqw = np.sum(IntegratedSqwg.reshape((-1, self.sqw.shape[0])),
                               axis=0)
        ene = np.arange(0, 3000, 1)
        spec = np.zeros(3000)
        for iw, s in enumerate(IntegratedSqw[1:]):
            dE = (self.E[iw+1] - self.E[0])
            #sigma = dE*0.02
            sigma = dE*0.04
            spec += s*np.exp(-(ene - dE)**2/sigma**2)
        if rank == 0:
            self.dataset = {}
            self.dataset['ene'] = ene
            self.dataset['spec'] = spec
            self.dataset['E'] = self.E
            self.dataset['integratedsqw'] = IntegratedSqw
        #    self.Plotter()
        if Iswrite and rank == 0:
            with open("./savedata_integrated.pkl", 'wb') as f:
                pickle.dump(self.dataset, f, 4)

