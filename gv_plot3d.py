#!/usr/bin/env python
# extract the reciplocal distribution of group velocities on a reciplocal plane.
# the output format is designed to be readable by VESTA, i.e., I gave up to make the 3D plot by myself.
# Written by Kazuyoshi TATSUMI
import numpy as np
import h5py
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import yaml

plt.style.use('classic')
plt.rcParams['font.family'] = 'Times New Roman'


def parse_mesh(filename):
    f = h5py.File(filename, 'r')
    gv = f["group_velocity"]
    qp = f["qpoint"][:, :]
    return gv, qp


def average_gv(filename, max_freq):
    print "filename", filename
    f = h5py.File(filename, 'r')
    if filename == "/home/kazu/bsi3n4_m/phono3py_doubled_112_fc2_334_sym_monk/kappa-m141416.bz.hdf5":
        print "chk", filename,"does not contain proper gv"
        gv, qpm = parse_mesh("/home/kazu/bsi3n4_m/phono3py_doubled_112_fc2_334_sym_monk/mesh_m141416.hdf5")
    else:
        gv = f["group_velocity"]

    gvdata = gv[:, :, 2]**2
    #gvdata = abs(gv[:, :, 2])
    #gvdata = gv[:, :, 2]**2 / (gv[:, :, 0]**2 + gv[:, :, 1]**2 + gv[:, :, 2]**2)
    #gvdata = (gv[:, :, 0])**2 + (gv[:, :, 1])**2
    omega = f["frequency"][:, :]
    mesh = f["mesh"]
    qp = f["qpoint"][:, :]
    qp[qp < 0] += 1
    dim = omega.shape
    _gvdata = np.zeros(dim)
    _gvdata[omega < max_freq] = gvdata[omega < max_freq]
    gvdata_ave = np.sum(_gvdata, axis=1)
    return(mesh, qp, gvdata_ave)


def average_kappa(filename, max_freq, temp):
    print filename
    f = h5py.File(filename, 'r')
    mesh = f["mesh"]
    print f.keys()
    kuc = f["kappa_unit_conversion"][()]
    qp = f["qpoint"][:, :]  # (gp, 3)
    print "qp",qp.shape
    if filename == "/home/kazu/bsi3n4_m/phono3py_doubled_112_fc2_334_sym_monk/kappa-m141416.bz.hdf5":
        print "chk", filename,"does not contain proper gv"
        gv, qpm = parse_mesh("/home/kazu/bsi3n4_m/phono3py_doubled_112_fc2_334_sym_monk/mesh_m141416.hdf5")
        qpm[qpm < 0] += 1
        print "qpm", qpm.shape
        print np.where((qp - qpm) > 0.01)
        print qpm[86,:], qp[86,:]
    else:
        gv = f["group_velocity"]
    #gv = f["group_velocity"] # (gp, band, 3)
    qp[qp < 0] += 1
    omega = f["frequency"][:, :] # (gp, band)
    gamma = f['gamma'][:] #KT, (temps, gp, band)
    cv = f['heat_capacity'][:] # (temps, gp, band)
    t = np.array(f['temperature'][:], dtype='double') # (temps)
    condition = abs(t - temp) < 0.001
    cv = np.squeeze(cv[condition, :, :])
    gamma = np.squeeze(gamma[condition, :, :])
    kappa = kuc * cv * gv[:, :, 2] * gv[:, :, 2] / (2 * gamma)
    kappa[0, 0:3] = 0
    #kappa = gv[:, :, 2]**2 
    print "kappa=",np.sum(kappa) / (mesh[0]*mesh[1]*mesh[2])
    dim = omega.shape
    _kappa = np.zeros(dim)
    _kappa[omega < max_freq] = kappa[omega < max_freq]
    kappa_ave = np.sum(_kappa, axis=1)
    print "kappa with omega < %s THz = %f" %  (max_freq, np.sum(kappa_ave) / (mesh[0]*mesh[1]*mesh[2]) )
    return(mesh, qp, kappa_ave)



def getXYZD(qx, qy, qz, gvdata_ave, n):
    xlin = np.linspace(min(qx), max(qx), n[0]*2)
    ylin = np.linspace(min(qy), max(qy), n[1]*2)
    zlin = np.linspace(min(qz), max(qz), n[2]*2)
    X, Y, Z = np.meshgrid(xlin, ylin, zlin)
    D = griddata((qx, qy, qz), gvdata_ave, (X, Y, Z), method='linear')

    D2 = np.transpose(D, (2, 1, 0))             # The output format is adjusted to the CHGCAR file.
    D3 = D2.reshape(n[0]*n[1]*n[2] / 5, 5)      # The transposing is nescessary for it. 
    return X, Y, Z, D3


def getXYZD2(q, kappadata, mesh):
    xlin = np.linspace(min(q[:, 0]), max(q[:, 0]), mesh[0])
    ylin = np.linspace(min(q[:, 1]), max(q[:, 1]), mesh[1])
    zlin = np.linspace(min(q[:, 2]), max(q[:, 2]), mesh[2])
    X, Y, Z = np.meshgrid(xlin, ylin, zlin)
    D = griddata((q[:, 0], q[:, 1], q[:, 2]), kappadata, (X, Y, Z), method='linear')

    return X, Y, Z, D


def getXYZ(x, y, gvdata_sum, n):
    xlin = np.linspace(min(x), max(x), n[0]*20)
    ylin = np.linspace(min(y), max(y), n[1]*20)
    X, Y = np.meshgrid(xlin, ylin)
    Z = griddata((x, y), gvdata_sum, (X, Y), method='linear')
    return X, Y, Z


def make_diff(file1, file2, dn, max_freq, rlat1, rlat2):
    mesh1, qp1, gvdata_ave1 = average_gv(file1, max_freq)
    mesh2, qp2, gvdata_ave2 = average_gv(file2, max_freq)
    X1, Y1, Z1, D1 = getXYZD(qp1[:, 0], qp1[:, 1], qp1[:, 2], gvdata_ave1, mesh1[:])
    X2, Y2, Z2, D2 = getXYZD(qp2[:, 0], qp2[:, 1], qp2[:, 2], gvdata_ave2, mesh2[:])
    diff = D2 * np.linalg.det(rlat2) - D1 * np.linalg.det(rlat1)
    return diff


def set_direction(mesh, qp, plane):
    if plane == "basal":
        meshx=mesh[0]
        meshy=mesh[1]
        meshz=mesh[2]
        qpx = qp[:, 0]
        qpy = qp[:, 1]
    elif plane == "prismatic":
        meshx=mesh[1]
        meshy=mesh[2]
        meshz=mesh[0]
        qpx = qp[:, 1]
        qpy = qp[:, 2]
    return meshx, meshy, meshz, qpx, qpy

def project_on_plane(mesh, qp, d, plane):
    k = 0
    meshx, meshy, meshz, qpx, qpy = set_direction(mesh, qp, plane)
    d_sum = np.zeros(((meshx+1)*(meshy+1)))
    x = np.zeros(((meshx+1)*(meshy+1)))
    y = np.zeros(((meshx+1)*(meshy+1)))
    for i in range(0, meshx):
        for j in range(0, meshy):
            x[k] = i * 1.0/meshx
            y[k] = j * 1.0/meshy
            condition = abs(qpx - x[k]) + abs(qpy - y[k]) < 0.001
            _d = d[condition]
            d_sum[k] = np.sum(_d)
            k += 1
    for i in range(0, meshx):
        x[k] = i * 1.0/meshx
        y[k] = 1.0
        condition = abs(qpx - x[k]) + abs(qpy - 0.000) < 0.001
        _d = d[condition]
        d_sum[k] = np.sum(_d)
        k += 1
    for i in range(0, meshy):
        x[k] = 1.0
        y[k] = i * 1.0/meshy
        condition = abs(qpx - 0.000) + abs(qpy - y[k]) < 0.001
        _d = d[condition]
        d_sum[k] = np.sum(_d)
        k += 1
    x[k] = 1.0
    y[k] = 1.0
    condition = abs(qpx - 0.000) + abs(qpy - 0.000) < 0.001
    _d = d[condition]
    d_sum[k] = np.sum(_d)
    return x, y, d_sum / meshz


def project_on_axis(mesh, qp, d, plane):
    k = 0
    meshx, meshy, meshz, qpx, qpy = set_direction(mesh, qp, plane)
    d_sum = np.zeros(((meshx+1)*(meshy+1)))
    x = np.zeros(((meshx+1)*(meshy+1)))
    y = np.zeros(((meshx+1)*(meshy+1)))
    for j in range(0, meshy):
        y[k] = j * 1.0/meshy
        condition = abs(qpy - y[k]) < 0.001
        _d = d[condition]
        d_sum[k] = np.sum(_d)
        k += 1
    y[k] = 1.0
    d_sum[k] = d_sum[0]
    return y, d_sum / (meshz * meshx)


def make_diffonplane(file1, file2, dn, max_freq, temp, rlat1, rlat2, plane):
    if dn == "group_velocity":
        mesh1, qp1, gvdata_ave1 = average_gv(file1, max_freq)
        print "gv_sum within omega < %s THz = %f" % (max_freq, np.sum(gvdata_ave1) / (mesh1[0]*mesh1[1]*mesh1[2]) * np.linalg.det(rlat1))
        x1, y1, z1 = project_on_plane(mesh1, qp1, gvdata_ave1*np.linalg.det(rlat1), plane)
        mesh2, qp2, gvdata_ave2 = average_gv(file2, max_freq)
        print "gv_sum within omega < %s THz = %f" % (max_freq, np.sum(gvdata_ave2) / (mesh2[0]*mesh2[1]*mesh2[2]) * np.linalg.det(rlat2))
        x2, y2, z2 = project_on_plane(mesh2, qp2, gvdata_ave2*np.linalg.det(rlat2), plane)
    if dn == "kappa":
        mesh1, qp1, kappadata_ave1 = average_kappa(file1, max_freq, temp)
        x1, y1, z1 = project_on_plane(mesh1, qp1, kappadata_ave1, plane)
        mesh2, qp2, kappadata_ave2 = average_kappa(file2, max_freq, temp)
        x2, y2, z2 = project_on_plane(mesh2, qp2, kappadata_ave2, plane)
        XX, YY, ZZ, DD = getXYZD2(qp1, (kappadata_ave2 - kappadata_ave1) / np.average(kappadata_ave2 - kappadata_ave1), mesh1)
        DDave = np.average(DD, axis=(0,1))
        DDave2 = np.tile(DDave, mesh1[0]*mesh1[1]).reshape(mesh1[0],mesh1[1],mesh1[2])
        plt.figure()
        plt.plot(DDave2[0, 0,:], ZZ[0,0,:])
        plt.plot(DDave, ZZ[0,0,:])

        SD = np.sum((DD - DDave2)**2, axis=2)
        loc = np.unravel_index(np.argmin(SD), SD.shape)
        plt.plot(DD[loc[0], loc[1], :], ZZ[0,0,:])
        

    if plane == "basal": 
        xy = np.zeros(((mesh2[0]+1)*(mesh2[1]+1), 2))
        xy[:, 0] = x2
        xy[:, 1] = y2
        carxy = np.matmul(xy, rlat2[0:2, 0:2])
    if plane == "prismatic": 
        xy = np.zeros(((mesh2[1]+1)*(mesh2[2]+1), 2))
        xy[:, 0] = x2
        xy[:, 1] = y2
        carxy = np.matmul(xy, rlat2[1:3, 1:3])
    X, Y, Z = getXYZ(carxy[:, 0], carxy[:, 1], (z2 - z1) / np.average(z2 - z1), mesh1[:])
    #Zave = np.average(Z[0:279,0:279], axis=1)
    #Zave = np.append(Zave, Zave[0])
    #plt.figure()
    #plt.plot(Zave, Y[:,0])
    #plt.plot(Z[:,0],Y[:,0])
    plt.figure()
    im = plt.contour(X, Y, Z, np.arange(0, 3, 0.1))
    plt.clabel(im, inline=1, fontsize=12, fmt='%3.1f')
    plt.title(dn)
    ax = plt.gca()
    ax.set_facecolor('k')


def parse_rlat(pf):
    with open(pf) as f:
        data = yaml.load(f)
        rlat = data["primitive_cell"]["reciprocal_lattice"]
    rlat = np.array(rlat)
    return rlat


def calc_modekappa(cv, gv, gamma, t, my_t):
    condition = abs(t - my_t) < 0.01
    print condition


def run():
    #cdir = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift"
    cdir = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/noiso"
    sdir = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/noiso"
    #sdir = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/noiso"
    ssdir = "/home/kazu/bsi3n4_m/phono3py_doubled_112_fc2_334_sym_monk"
    #ssdir = "/home/kazu/bsi3n4_m/phonopy_doubled_334"
    #ssdir = "/home/kazu/bsi3n4_m/phonopy_doubled_334_with_alphalat"
    max_freq = 45
    temp = 300
    #pc = cdir + "/phonopy.yaml"
    pc = cdir + "/phono3py.yaml"
    #ps = sdir + "/phonopy.yaml"
    ps = sdir + "/phono3py.yaml"
    #pss = ssdir + "/phonopy.yaml"
    pss = ssdir + "/phono3py.yaml"
    crlat = parse_rlat(pc)
    srlat = parse_rlat(ps)
    ssrlat = parse_rlat(pss)
    #c = cdir + "/mesh_252535.hdf5"
    #c = cdir + "/mesh.hdf5.org"
    c = cdir + "/kappa-m141416.bz.hdf5"
    #s = sdir + "/mesh.hdf5.gonze"
    s = sdir + "/kappa-m141432.bz.hdf5"
    #s = sdir + "/gpjob_m101026_with_recent_version/kappa-m101026.bz.hdf5"    #kappa file with data over all mesh points, generated by expand_ibz_2_solid.py
    ss = ssdir + "/kappa-m141416.bz.hdf5"    #kappa file with data over all mesh points, generated by expand_ibz_2_solid.py
    #ss = ssdir + "/mesh_m141416.hdf5"    #kappa file with data over all mesh points, generated by expand_ibz_2_solid.py
    #ss = ssdir + "/mesh_252535.hdf5"
    #ss = ssdir + "/mesh_m101014.hdf5"
    ##For data extraction for each phase
    #mesh, qp, gvdata_ave = average_gv(c, max_freq)
    #X, Y, Z, D3 = getXYZD(qp[:, 0], qp[:, 1], qp[:, 2], gvdata_ave, mesh[:])
    #np.savetxt("array.txt", D3 * np.linalg.det(srlat), fmt="%17.11e")  #
    ##Seek the difference between alpha and beta_doubled_cell
    #diff = make_diff(c, ss, "group_velocity", max_freq, crlat, ssrlat)
    #plt.figure()
    #make_diffonplane(c, ss, "kappa", max_freq, temp, crlat, ssrlat, "prismatic")
    make_diffonplane(c, ss, "kappa", max_freq, temp, crlat, ssrlat, "prismatic")
    #plt.figure()
    make_diffonplane(c, ss, "group_velocity", max_freq, temp, crlat, ssrlat, "prismatic")
    #make_diffonplane(c, s, "kappa", max_freq, temp, crlat, srlat, "basal")
    #np.savetxt("array.txt", diff, fmt="%17.11e")  #
    plt.show()



run()
