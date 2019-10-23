#!/usr/bin/env python

import numpy as np
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def get_phonon_data(pre_dir):
    f = h5py.File(pre_dir + "mesh.hdf5", 'r')
    qpoints = f['qpoint'][:]
    freqs = f['frequency'][:]
    eigvecs = f['eigenvector'][:]
    gvs = f['group_velocity'][:]

    return qpoints, freqs, eigvecs, gvs


def get_kappa_data(pre_dir, filename):
    f = h5py.File(pre_dir + "noiso/" + filename, 'r')
    gamma = f['gamma'][:]
    qpoints = f['qpoint'][:]
    return gamma, qpoints


def check_qpoints(q_ir, q_all, grid_mapping):
    grid_unique = np.unique(grid_mapping)
    q_unique = q_all[grid_unique]
    assert len(q_ir) == len(q_unique)

    for q1, q2 in zip(q_ir, q_unique):
        print(q1, q2)

    for q in q_unique:
        diff = q_ir - q
        diff -= np.rint(diff)
        match = np.where((diff ** 2).sum(axis=1) < 1e-5)[0]
        assert len(match) == 1, "%s" % diff

    for q1, q2 in zip(q_ir, q_unique):
        diff = q1 - q2
        diff -= np.rint(diff)
        assert (diff ** 2).sum() < 1e-5


def plot(lifetime, freqs, eigvecs):
    fig, ax = plt.subplots()
    natom = freqs.shape[1] // 3
    idx = [3 * i for i in range(natom)]
    e_z = eigvecs[:, idx, :]
    e_z2 = (np.abs(e_z) ** 2).sum(axis=1)

    df = np.random.rand(*(freqs.shape)) / 20
    dl = np.random.rand(*(freqs.shape))

    sc = ax.scatter(freqs + df, lifetime + dl, s=4, c=e_z2, cmap='PiYG')
    ax.set_xlim(0, 35)
    fig.colorbar(sc)
    x = 0.1*np.arange(0, 51)
    y = (-55*x + 515) / (2 * np.pi)
    ax.plot(x, y)


def plotgv(gvs, freqs):
    fig, ax = plt.subplots(2, 2)

    df = np.random.rand(*(freqs.shape)) / 20
    dl = np.random.rand(*(freqs.shape))
    gv_mag = (np.abs(gvs) ** 2).sum(axis=2)
    gv_z = np.abs(gvs[:, :, 2]) ** 2

    sc = ax[0, 0].scatter(freqs + df, gv_z + dl, s=1)
    ax[0, 0].set_xlim(0, 15)
    ax[0, 0].set_ylim(0, 17000)
    sc = ax[0, 1].scatter(freqs + df, gv_mag + dl, s=1)
    ax[0, 1].set_xlim(0, 15)
    ax[0, 1].set_ylim(0, 17000)


def plotq(lifetime, freqs, qpoints_all, eigvecs):
    fig = plt.figure(figsize=(16, 6))
    for mindx in range(0, 3):
        arr, arr2, earr, earr2 = get_mode(lifetime, freqs, qpoints_all, eigvecs, mindx)

        natom = earr.shape[1] // 3
        idx = [3 * i for i in range(natom)]
        e_z = earr[:, idx]
        e2_z = earr2[:, idx]
        e_z2 = (np.abs(e_z) ** 2).sum(axis=1)
        e2_z2 = (np.abs(e2_z) ** 2).sum(axis=1)

        ax = fig.add_subplot(2, 3, mindx + 1, projection='3d')
        ax2 = fig.add_subplot(2, 3, mindx + 4, projection='3d')
        #df = np.random.uniform(-1, 1, (arr[:, 0].shape)) / 80
        #dl = np.random.uniform(-1, 1, (arr[:, 0].shape)) / 80
        #dv = np.random.uniform(-1, 1, (arr[:, 0].shape)) / 80
        #df2 = np.random.uniform(-1, 1, (arr2[:, 0].shape)) / 80
        #dl2 = np.random.uniform(-1, 1, (arr2[:, 0].shape)) / 80
        #dv2 = np.random.uniform(-1, 1, (arr2[:, 0].shape)) / 80

        #ax.scatter3D(arr[:, 0] + df, arr[:, 1] + dl, arr[:, 2] + dv, c=e_z2, cmap='PiYG', s=4)
        #ax2.scatter3D(arr2[:, 0] + df2, arr2[:, 1] + dl2, arr2[:, 2] + dv2, c=e2_z2, cmap='PiYG',  s=4)
        #ax.scatter3D(arr[:, 0], arr[:, 1], arr[:, 2], c=e_z2, cmap='PiYG', s=4)
        p=ax.scatter3D(arr2[:, 0], arr2[:, 1], arr2[:, 2], c=e2_z2, cmap='PiYG', s=4)
        p2 = ax2.scatter3D(arr[:, 0], arr[:, 1], arr[:, 2], c=e_z2, cmap='PiYG', s=4)
        ax.view_init(azim=0, elev=0)
        ax2.view_init(azim=0, elev=0)
        ax.set_xlabel('qx')
        ax.set_ylabel('qy')
        ax.set_zlabel('qz')
        ax.set_zlim(-0.50, 0.50)
        ax2.set_zlim(-0.50, 0.50)
        ax.set_xlim(-0.6, 0.6)
        ax2.set_xlim(-0.6, 0.6)
        ax.set_ylim(-0.6, 0.6)
        ax2.set_ylim(-0.6, 0.6)
        fig.colorbar(p2, ax=ax2)
        fig.colorbar(p, ax=ax)


def plot_alpha():
    pre_dir = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/"
    qpoints_all, freqs, eigvecs, gvs = get_phonon_data(pre_dir)
    f = h5py.File(pre_dir + "grid_address-m101014.hdf5", 'r')
    grid_mapping = f['grid_mapping_table'][:]
    gamma_ir, qpoints_ir = get_kappa_data(pre_dir,
                                          "kappa-m101014.noiso.hdf5")
    check_qpoints(qpoints_ir, qpoints_all, grid_mapping)
    gamma_at_300K_ir = gamma_ir[20]
    g_ir = np.where(gamma_at_300K_ir > 0, gamma_at_300K_ir, -1)
    lifetime_ir = np.where(g_ir > 0, 1.0 / (2 * 2 * np.pi * g_ir), 0)
    lifetime = get_lifetime_all(lifetime_ir, grid_mapping)
    #plot(lifetime, freqs, eigvecs)
    plotgv(gvs, freqs)


def get_lifetime_all(lifetime_ir, grid_mapping):
    ngp = len(grid_mapping)
    lifetime = np.zeros(((ngp, ) + lifetime_ir.shape[1:]),
                        dtype=lifetime_ir.dtype)
    unique_gp = np.unique(grid_mapping)
    gp_map = {gp: i for i, gp in enumerate(unique_gp)}

    for i, gp in enumerate(grid_mapping):
        lifetime[i] = lifetime_ir[gp_map[gp]]
    return lifetime


def get_mode(lifetime, freqs, qpoints_all, eigvecs, mindx):
    diff_lifetime = lifetime - (-55*freqs + 515) / (2 * np.pi)
    match = np.where(diff_lifetime > 0)
    match_qindx = match[0]
    match_pindx = match[1]
    unmatch = np.where(diff_lifetime <= 0)
    unmatch_qindx = unmatch[0]
    unmatch_pindx = unmatch[1]
    #idx = [3 * i for i in range(natom)]
    #e_z = eigvecs[:, idx, :]
    #e_z2 = (np.abs(e_z) ** 2).sum(axis=1)
    
    #match = np.where(e_z2 < 0.1)
    #match_qindx = match[0]
    #match_pindx = match[1]
    #unmatch = np.where(e_z2 >= 0.1)
    #unmatch_qindx = unmatch[0]
    #unmatch_pindx = unmatch[1]

    arr = np.empty((0, 3), float)
    arr2 = np.empty((0, 3), float)
    earr = np.empty((0, 42), float)
    earr2 = np.empty((0, 42), float)
    for i, mq in enumerate(match_qindx):
        if freqs[match_qindx[i], match_pindx[i]] < 15 and match_pindx[i] == mindx:
            #print qpoints_all[match_qindx[i]], match_pindx[i]
            arr = np.concatenate((arr, [qpoints_all[match_qindx[i]]]), axis=0)
            earr = np.concatenate((earr, [eigvecs[match_qindx[i], :, match_pindx[i]]]), axis=0)
    for i, mq in enumerate(unmatch_qindx):
        if freqs[unmatch_qindx[i], unmatch_pindx[i]] < 15 and unmatch_pindx[i] == mindx:
            arr2 = np.concatenate((arr2, [qpoints_all[unmatch_qindx[i]]]), axis=0)
            earr2 = np.concatenate((earr2, [eigvecs[unmatch_qindx[i], :, unmatch_pindx[i]]]), axis=0)
    return arr, arr2, earr, earr2


def plot_beta():
    pre_dir = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/"
    qpoints_all, freqs, eigvecs, gvs = get_phonon_data(pre_dir)
    f = h5py.File(pre_dir + "grid_address-m101026.hdf5", 'r')
    grid_mapping = f['grid_mapping_table'][:]
    gamma_ir, qpoints_ir = get_kappa_data(pre_dir,
                                          "kappa-m101026.noiso.hdf5")
    check_qpoints(qpoints_ir, qpoints_all, grid_mapping)
    gamma_at_300K_ir = gamma_ir[20]
    g_ir = np.where(gamma_at_300K_ir > 0, gamma_at_300K_ir, -1)
    lifetime_ir = np.where(g_ir > 0, 1.0 / (2 * 2 * np.pi * g_ir), 0)
    lifetime = get_lifetime_all(lifetime_ir, grid_mapping)
    #plot(lifetime, freqs)
    plotgv(gvs, freqs)


def main():
    plot_alpha()
    plot_beta()
    plt.show()


if __name__ == '__main__':
    main()
