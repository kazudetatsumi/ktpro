#!/usr/bin/env python
# This script compares the densities estimated by the mem and sm methods.
# For verification of sparse modeling on the nuclear reconstruction proposed by 
# professor H. Tanaka. 
# Kazuyoshi TATSUMI 2022/02/10
import sys
sys.path.append("/home/kazu/ktpro")
import create_den as cd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import colors
import numpy as np


def set_axis(ax, cntr):
    ax.set_xlabel('z position (Angs)')
    ax.set_ylabel('x position (Angs)')
    ax.clabel(cntr, fmt='%1.0f', fontsize=6)
    ax.set_aspect('equal')


def get_integrated_density(xi, yi, zi, density, rcs):
    ctr = np.unravel_index(np.argmax(density, axis=None), density.shape)
    cpos = np.array([xi[ctr], yi[ctr], zi[ctr]])
    dist = ((xi - cpos[0])**2 + (yi - cpos[1])**2 + (zi - cpos[2])**2)**0.5
    rdist = np.zeros_like(rcs)
    for ridx,  rc in enumerate(rcs):
        rdist[ridx] = np.sum(density[dist < rc])
    return(rdist)


def check_density():
    mesh = 128
    latconst = 10.0
    atomd = 1.1
    atomw = 0.2
    outfile = 'input.den'
    protfile = 'prot.den'
    proj = cd.create_den(outfile, protfile, mesh, latconst, atomd, atomw)
    #proj.readden('dnsity.den')
    input_density = proj.readgrd('dnsity.grd')
    #print(np.unravel_index(np.argmax(input_density, axis=None), input_density.shape))
    sm_density = proj.readgrd('result/SM/out2_sm.grd')
    mem_density = proj.readgrd('result/MEM/out2_mem.grd')
    fig = plt.figure(figsize=(14, 24))
    gs = gridspec.GridSpec(nrows=5, ncols=3)
    gs1 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs[0, :])
    gs2 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs[1, :])
    gs21 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs[2, :])
    gs22 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs[3, :])
    gs3 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=3, subplot_spec=gs[4, :])
    #print(np.unravel_index(np.argmax(sm_density, axis=None), sm_density.shape))
    #print(np.unravel_index(np.argmax(mem_density, axis=None), mem_density.shape))

    # density on the prinicple axis
    ax = fig.add_subplot(gs1[0, :])
    plt.plot(np.linspace(0, 10, 128), input_density[64, 64, :], label='true')
    plt.plot(np.linspace(0, 10, 128), sm_density[64, 64, :], label='sm')
    plt.plot(np.linspace(0, 10, 128), mem_density[64, 64, :], label='mem')
    plt.xlabel('z position (Angs)')
    plt.ylabel('scattering length density (arb. units)')
    plt.xlim(2, 8)
    plt.legend()
    #plt.ylim(0, 100)
    u = np.linspace(2.5, 7.5, 64)
    v = np.linspace(2.5, 7.5, 64)
    X, Z = np.meshgrid(u, v)
    uf = np.linspace(0.0, 10.0, 128)
    vf = np.linspace(0.0, 10.0, 128)
    Xf, Zf = np.meshgrid(uf, vf)
    levels = np.power(np.arange(4), 6)

    # density on a plane containing the atom centers
    ax = fig.add_subplot(gs2[0, 0])
    plt.title('true')
    cntr = ax.contour(X, Z, input_density[32:96, 64, 32:96], levels, colors='black')
    set_axis(ax, cntr)
    ax = fig.add_subplot(gs2[0, 1])
    plt.title('sm')
    cntr = ax.contour(X, Z, sm_density[32:96, 64, 32:96], levels, colors='black')
    set_axis(ax, cntr)
    ax = fig.add_subplot(gs2[0, 2])
    plt.title('mem')
    cntr = ax.contour(X, Z, mem_density[32:96, 64, 32:96], levels, colors='black')
    set_axis(ax, cntr)

    # density on a plane containing the atom centers, gray scale is set for the
    # positions apart from the atom centers.
    ax = fig.add_subplot(gs21[0, 0])
    vmin = np.min(mem_density)/50000.0
    vmax = np.max(mem_density)/50000.0
    plt.title('true')
    c = ax.pcolor(Xf, Zf, input_density[:, 64, :], vmin=vmin, vmax=vmax,
                  shading='auto', cmap='gray_r')
    ax.set_aspect('equal')
    fig.colorbar(c, ax=ax)
    ax = fig.add_subplot(gs21[0, 1])
    c = ax.pcolor(Xf, Zf, sm_density[:, 64, :], vmin=vmin, vmax=vmax,
                  shading='auto', cmap='gray_r')
    ax.set_aspect('equal')
    fig.colorbar(c, ax=ax)
    plt.title('sm')
    ax = fig.add_subplot(gs21[0, 2])
    c = ax.pcolor(Xf, Zf, mem_density[:, 64, :], vmin=vmin, vmax=vmax,
                  shading='auto', cmap='gray_r')
    plt.title('mem')
    ax.set_aspect('equal')
    fig.colorbar(c, ax=ax)

    # We subtracted the true density from the extimated density and the
    # resultand difference densities are shown on a plane containing the atom
    # centers, color scale is set for the positions near the atom centers.
    ax = fig.add_subplot(gs22[0, 1])
    plt.title('sm-true')
    vmin = np.min(mem_density[32:96, 64, 32:96] -
                  input_density[32:96, 64, 32:96])
    vmax = np.max(mem_density[32:96, 64, 32:96] -
                  input_density[32:96, 64, 32:96])
    norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

    c = ax.pcolor(X, Z, sm_density[32:96, 64, 32:96] -
                  input_density[32:96, 64, 32:96], norm=norm, shading='auto',
                  cmap='bwr')
    ax.set_aspect('equal')
    fig.colorbar(c, ax=ax)
    ax = fig.add_subplot(gs22[0, 2])
    plt.title('mem-true')
    c = ax.pcolor(X, Z, mem_density[32:96, 64, 32:96] -
                  input_density[32:96, 64, 32:96], norm=norm, shading='auto',
                  cmap='bwr')
    ax.set_aspect('equal')
    fig.colorbar(c, ax=ax)

    # cumulative sum of scattering lengths
    xi, yi, zi = np.mgrid[0.0:10.0:128j, 0.0:10.0:128j, 0.0:10.0:128j]
    rc = np.linspace(0, 8.0, 160)
    input_rdist = get_integrated_density(xi, yi, zi, input_density, rc)
    sm_rdist = get_integrated_density(xi, yi, zi, sm_density, rc)
    mem_rdist = get_integrated_density(xi, yi, zi, mem_density, rc)
    #print(input_rdist.shape)
    #print(sm_rdist.shape)
    #print(mem_rdist.shape)
    ax = fig.add_subplot(gs3[0, :])
    plt.plot(rc, input_rdist, label='true')
    plt.plot(rc, sm_rdist, label='sm')
    plt.plot(rc, mem_rdist, label='mem')
    plt.xlabel('distance from the maximum density position (Angs)')
    plt.ylabel('cumulative sum of scattering lengths \n in the unit cell (arb. units)')
    plt.legend()
    plt.savefig("check_density.png")
    plt.show()
    #proj.get_all_coord()
    #proj.get_density()
    #proj.output()


check_density()
