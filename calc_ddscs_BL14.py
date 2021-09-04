#!/usr/bin/env python
# This script is a test script to calculate phonon mode-specific scattering cross-sections from a phonopy ouptupt hdf5 file.
# The phonopy output hdf5 file should contain "eigenvector", "frequency" and "qpoint".
# Atomic masses and coherent scattering lengths should be hard-coded in this script.
import h5py
import numpy as np
import math

temperature = 15                      # [K]
PlanckConstant = 4.13566733e-15       # [eV s]
THzTomev = PlanckConstant * 1e15      # [meV]
kb = 8.617333262145e-2                # [meV K-1]
Ecut = 0.1                           # [meV]
u2 = 0.0018364                          # [angs2]
Ein = 42.05                           # [meV]
rlat = 0.275360831107482          # [Angs-1]


def getphonopydatafromqpoint(fn, mesh):
    f = h5py.File(fn, 'r')
    nmodes = f["frequency"].shape[1]                                    # f["frequency"]:[number of qpoints, number of modes], ordered as z indx increasing faster
    omega = np.reshape(f["frequency"][:],
                       [mesh[0], mesh[1], mesh[2], nmodes], order='C')
    omega[omega < 0] = 0.0
    eigvec = np.reshape(f["eigenvector"][:],
            [mesh[0], mesh[1], mesh[2], nmodes, nmodes], order='C') # f["eigenvector"]:[number of qpoints, number of elements in eigenveoctor, number of modes]
    qvec = np.reshape(f["qpoint"], [mesh[0], mesh[1], mesh[2], 3], order='C')
    print(np.max(qvec))
    print(np.min(qvec))
    return(omega, eigvec, qvec)


def bose_distribution_func(omega, temperature):
    nb = 1.0 / (np.exp(omega/(temperature*kb)) - 1.0)
    return(nb)


def sumofxd_func(eigvec, qvec):
    return np.sum(qvec[:, :, :, :, np.newaxis]*eigvec, axis=3)


def save_h5py(ddscs, outfile):
    with h5py.File(outfile, 'w') as hf:
        hf.create_dataset('ddscs', data=ddscs)


def run():
    dqx = 0.0125
    dqy = 0.025
    dqz = 0.050
    lb_qx = -0.675
    ub_qx = 3.075 + dqx
    lb_qy = -0.925
    ub_qy = 4.375 + dqy
    lb_qz = -0.80
    ub_qz = 0.55 + dqz
    xlin = np.arange(lb_qx, ub_qx, dqx)
    nx = xlin.shape[0]
    ylin = np.arange(lb_qy, ub_qy, dqy)
    ny = ylin.shape[0]
    zlin = np.arange(lb_qz, ub_qz, dqz)
    nz = zlin.shape[0]
    mesh = np.array([nx, ny, nz])

    qpointsfile = "/home/kazu/WORK/vasp-phonopy/cu_BL14/test_kfki_debye/qpoints.hdf5"
    outfile = "/home/kazu/WORK/vasp-phonopy/cu_BL14/test_kfki_debye/ddscs.hdf5"
    omega, eigvec, qvec = getphonopydatafromqpoint(qpointsfile, mesh)
    print(omega.shape)
    print(eigvec.shape)
    print(qvec.shape)
    omega = omega*THzTomev
    nb = bose_distribution_func(omega, temperature)
    sumofxd = sumofxd_func(eigvec, qvec)
    kfki = np.sqrt((Ein - omega)/Ein)
    _qvec = qvec * rlat * 2.0 * math.pi
    dwfac = np.zeros(omega.shape)
    for imode in range(0, dwfac.shape[3]):
        dwfac[:, :, :, imode] = np.exp(-np.sum(_qvec*_qvec, axis=3)*u2)
    ddscs = kfki * dwfac * np.abs(sumofxd)**2 * (nb + 1.0) / omega
    ddscs[omega < Ecut] = 0.
    save_h5py(ddscs, outfile)


run()
