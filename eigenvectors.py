#!/usr/bin/env python

import numpy as np
import h5py
import matplotlib.pyplot as plt


def get_phonon_data(pre_dir):
    f = h5py.File(pre_dir + "mesh.hdf5", 'r')
    qpoints = f['qpoint'][:]
    freqs = f['frequency'][:]
    eigvecs = f['eigenvector'][:]

    return qpoints, freqs, eigvecs


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
    ax.set_xlim(0, 5)
    fig.colorbar(sc)
    x = 0.1*np.arange(0, 51)
    y = (-55*x + 515) / (2 * np.pi)
    ax.plot(x, y)

    #print e_z2.shape
    #diff_lifetime = lifetime - (-55*freqs + 515) / (2 * np.pi)
    #match = np.where(diff_lifetime > 0)
    #match_qindx = match[0]
    #match_pindx = match[1]
    #counter = 0 
    #for i, mq in enumerate(match_qindx):
    #    if freqs[match_qindx[i], match_pindx[i]] < 5:
    #        print qpoints_all[match_qindx[i]], match_pindx[i]
    #        counter += 1
    #print counter


def plot_alpha():
    pre_dir = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/"
    qpoints_all, freqs, eigvecs = get_phonon_data(pre_dir)
    f = h5py.File(pre_dir + "grid_address-m101014.hdf5", 'r')
    grid_mapping = f['grid_mapping_table'][:]
    gamma_ir, qpoints_ir = get_kappa_data(pre_dir,
                                          "kappa-m101014.noiso.hdf5")
    check_qpoints(qpoints_ir, qpoints_all, grid_mapping)
    gamma_at_300K_ir = gamma_ir[20]
    g_ir = np.where(gamma_at_300K_ir > 0, gamma_at_300K_ir, -1)
    lifetime_ir = np.where(g_ir > 0, 1.0 / (2 * 2 * np.pi * g_ir), 0)
    lifetime = get_lifetime_all(lifetime_ir, grid_mapping)
    plot(lifetime, freqs, eigvecs)


def get_lifetime_all(lifetime_ir, grid_mapping):
    ngp = len(grid_mapping)
    lifetime = np.zeros(((ngp, ) + lifetime_ir.shape[1:]),
                        dtype=lifetime_ir.dtype)
    unique_gp = np.unique(grid_mapping)
    gp_map = {gp: i for i, gp in enumerate(unique_gp)}

    for i, gp in enumerate(grid_mapping):
        lifetime[i] = lifetime_ir[gp_map[gp]]
    return lifetime


def plot_beta():
    pre_dir = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/"
    qpoints_all, freqs, eigvecs = get_phonon_data(pre_dir)
    f = h5py.File(pre_dir + "grid_address-m101026.hdf5", 'r')
    grid_mapping = f['grid_mapping_table'][:]
    gamma_ir, qpoints_ir = get_kappa_data(pre_dir,
                                          "kappa-m101026.noiso.hdf5")
    check_qpoints(qpoints_ir, qpoints_all, grid_mapping)
    gamma_at_300K_ir = gamma_ir[20]
    g_ir = np.where(gamma_at_300K_ir > 0, gamma_at_300K_ir, -1)
    lifetime_ir = np.where(g_ir > 0, 1.0 / (2 * 2 * np.pi * g_ir), 0)
    lifetime = get_lifetime_all(lifetime_ir, grid_mapping)
    plot(lifetime, freqs, eigvecs)

    print eigvecs.shape
    diff_lifetime = lifetime - (-55*freqs + 515) / (2 * np.pi)
    match = np.where(diff_lifetime > 0)
    match_qindx = match[0]
    match_pindx = match[1]
    counter = 0 
    for i, mq in enumerate(match_qindx):
        if freqs[match_qindx[i], match_pindx[i]] < 5:
            print qpoints_all[match_qindx[i]], match_pindx[i]
            counter += 1
    print counter
    

def main():
    plot_beta()
    plt.show()
    


if __name__ == '__main__':
    main()
