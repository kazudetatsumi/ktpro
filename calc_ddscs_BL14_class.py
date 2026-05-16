#!/usr/bin/env python
# This script is a test script to calculate phonon mode-specific scattering
# cross-sections from a phonopy ouptupt hdf5 file.
# The phonopy output hdf5 file should contain "eigenvector", "frequency" and
# "qpoint".
# Atomic masses and coherent scattering lengths should be hard-coded in this
# script.
import h5py
import numpy as np
import math


class CALC_ddscs:
    def __init__(self, qpointsfile, outfile,
                 temperature=15, Ecut=0.1, u2=0.0018364, Ein=42.05,
                 rlat=0.275360831107482, dqx=0.0125, dqy=0.025, dqz=0.050,
                 lb_qx=-0.675, ub_qx=3.075, lb_qy=-0.925, ub_qy=4.375,
                 lb_qz=-0.80, ub_qz=0.55):
        self.temperature = temperature             # [K]
        self.PlanckConstant = 4.13566733e-15       # [eV s]
        self.THzTomev = self.PlanckConstant * 1e15      # [meV]
        self.kb = 8.617333262145e-2                # [meV K-1]
        self.Ecut = 0.1                            # [meV]
        self.u2 = u2                               # [angs2]
        self.Ein = Ein                             # [meV]
        self.rlat = rlat                           # [Angs-1]
        ub_qx += dqx
        ub_qy += dqy
        ub_qz += dqz
        xlin = np.arange(lb_qx, ub_qx, dqx)
        nx = xlin.shape[0]
        ylin = np.arange(lb_qy, ub_qy, dqy)
        ny = ylin.shape[0]
        zlin = np.arange(lb_qz, ub_qz, dqz)
        nz = zlin.shape[0]
        self.mesh = np.array([nx, ny, nz])
        self.qpointsfile = qpointsfile
        self.outfile = outfile
        print("inverse temperature:", 1/(self.kb*self.temperature), "[meV^-1]")

    def getphonopydatafromqpoint(self):
        f = h5py.File(self.qpointsfile, 'r')
        # f["frequency"]:[number of qpoints, number of modes], ordered as z
        # indx increasing faster
        nmodes = f["frequency"].shape[1]
        self.omega = np.reshape(f["frequency"][:], [self.mesh[0], self.mesh[1],
                                self.mesh[2], nmodes], order='C')
        self.omega[self.omega < 0] = 0.0
        self.omega = self.omega*self.THzTomev
        # f["eigenvector"]:[number of qpoints, number of elements in
        # eigenveoctor, number of modes]
        self.eigvec = np.reshape(f["eigenvector"][:], [self.mesh[0],
                                 self.mesh[1], self.mesh[2], nmodes, nmodes],
                                 order='C')
        self.qvec = np.reshape(f["qpoint"], [self.mesh[0], self.mesh[1],
                               self.mesh[2], 3], order='C')
        print(np.max(self.qvec))
        print(np.min(self.qvec))

    def bose_distribution_func(self):
        self.nb = 1.0 / (np.exp(self.omega/(self.temperature*self.kb)) - 1.0)

    def sumofxd_func(self):
        self.sumofxd = np.sum(self.qvec[:, :, :, :, np.newaxis]*self.eigvec,
                              axis=3)

    def save_h5py(self):
        with h5py.File(self.outfile, 'w') as hf:
            hf.create_dataset('ddscs', data=self.ddscs)

    def honegumi(self):
        self.getphonopydatafromqpoint()
        print(self.omega.shape)
        print(self.eigvec.shape)
        print(self.qvec.shape)
        self.bose_distribution_func()
        self.sumofxd_func()
        self.kfki = np.sqrt((self.Ein - self.omega)/self.Ein)
        _qvec = self.qvec * self.rlat * 2.0 * math.pi
        dwfac = np.exp(-np.sum(_qvec*_qvec, axis=3)*self.u2)[:, :, :, np.newaxis]
        self.ddscs = self.kfki * dwfac * np.abs(self.sumofxd)**2 *\
            (self.nb + 1.0) / self.omega
        self.ddscs[self.omega < self.Ecut] = 0.
        #self.save_h5py()

    def get_banddisp(self):
        self.getphonopydatafromqpoint()
        print(self.omega.shape)
        print(self.qvec.shape)


def samplerun():
    dqx = 0.0125
    dqy = 0.025
    dqz = 0.050
    lb_qx = -0.675
    ub_qx = 3.075
    lb_qy = -0.925
    ub_qy = 4.375
    lb_qz = -0.80
    ub_qz = 0.55
    qpointsfile = "/home/kazu/WORK/vasp-phonopy/cu_BL14/test_kfki_debye/" +\
                  "qpoints.hdf5"
    outfile = "/home/kazu/WORK/vasp-phonopy/cu_BL14/test_kfki_debye/ddscs.hdf5"
    proj = CALC_ddscs(qpointsfile, outfile, dqx=dqx, dqy=dqy, dqz=dqz, lb_qx=lb_qx, ub_qx=ub_qx,
                      lb_qy=lb_qy, ub_qy=ub_qy, lb_qz=lb_qz, ub_qz=ub_qz)
    proj.honegumi()
    print(proj.ddscs)


#samplerun()
