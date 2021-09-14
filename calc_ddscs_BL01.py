#!/usr/bin/env python
# This script is a test script to calculate phonon mode-specific scattering cross-sections from a phonopy ouptupt hdf5 file.
# The phonopy output hdf5 file should contain "eigenvector", "frequency" and "qpoint".
# Atomic masses and coherent scattering lengths should be hard-coded in this script.
import h5py
import numpy as np
import math

temperature = 292                     # [K]
PlanckConstant = 4.13566733e-15       # [eV s]
THzTomev = PlanckConstant * 1e15      # [meV]
kb = 8.617333262145e-2                # [meV K-1]
Ecut = 0.1                            # [meV]
u2 = 0.00623                          # [angs2]
Ein = 50.0                           # [meV]
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
    return(omega, eigvec, qvec)


def bose_distribution_func(omega, temperature):
    nb = 1.0 / (np.exp(omega/(temperature*kb)) - 1.0)
    return(nb)


def sumofxd_func(eigvec, qvec, omega):
    # at present, assuming a Bravais Crystal, where the atom position is only (0,0,0).
    xd = np.zeros((eigvec.shape[0], eigvec.shape[1], eigvec.shape[2], eigvec.shape[4]), dtype=np.complex)
    for imode in range(0, eigvec.shape[4]):
        xd[:, :, :, imode] = (qvec[:, :, :, 0]*eigvec[:, :, :, 0, imode] +
                              qvec[:, :, :, 1]*eigvec[:, :, :, 1, imode] +
                              qvec[:, :, :, 2]*eigvec[:, :, :, 2, imode])
    return(xd)


def save_h5py(ddscs, outfile):
    with h5py.File(outfile, 'w') as hf:
        hf.create_dataset('ddscs', data=ddscs)


def run():
    dqx = 0.025
    dqy = 0.025
    dqz = 0.025
    lb_qx = -1.65
    ub_qx = 4.10 
    lb_qy = -2.1
    ub_qy = 2.8 
    lb_qz = -0.85
    ub_qz = 0.9 
    xlin = np.arange(lb_qx, ub_qx, dqx)
    nx = xlin.shape[0]
    ylin = np.arange(lb_qy, ub_qy, dqy)
    ny = ylin.shape[0]
    zlin = np.arange(lb_qz, ub_qz, dqz)
    nz = zlin.shape[0]
    mesh = np.array([nx, ny, nz])
    print(mesh)

    qpointsfile = "/home/kazu/WORK/vasp-phonopy/cu_BL01/test_kfki_debye_rt/qpoints.hdf5"
    outfile = "/home/kazu/WORK/vasp-phonopy/cu_BL01/test_kfki_debye_rt/ddscs.hdf5"
    omega, eigvec, qvec = getphonopydatafromqpoint(qpointsfile, mesh)
    print(omega.shape)
    print(eigvec.shape)
    print(qvec.shape)
    omega = omega*THzTomev
    nb = bose_distribution_func(omega, temperature)
    sumofxd = sumofxd_func(eigvec, qvec, omega)
    kfki = np.sqrt((Ein - omega)/Ein)
    _qvec = qvec * rlat * 2.0 * math.pi
    dwfac = np.zeros(omega.shape)
    for imode in range(0, dwfac.shape[3]):
        dwfac[:, :, :, imode] = np.exp(-np.sum(_qvec*_qvec, axis=3)*u2)
    ddscs = kfki * dwfac * np.abs(sumofxd)**2 * (nb + 1.0) / omega
    ddscs[omega < Ecut] = 0.
    print(np.sum(ddscs))
    save_h5py(ddscs, outfile)


run()



