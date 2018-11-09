#!/usr/bin/env python

from matplotlib.colors import LinearSegmentedColormap
import glob
import h5py
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.cm as cm

plt.style.use('classic')
plt.rcParams['font.family'] = 'Times New Roman'


def generate_cmap(colors):
    values = range(len(colors))

    vmax = np.ceil(max(values))
    color_list = []
    for v, c in zip(values, colors):
        color_list.append((v / vmax, c))
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)


def contour(i, freq, x, y, xt, yt, pn):
    ny = np.ceil((max(y)-min(y))*1200)
    nx = np.ceil((max(x)-min(x))*1200)
    print ny, nx
    xlin = np.linspace(min(x), max(x), nx)
    ylin = np.linspace(min(y), max(y), ny)
    X, Y = np.meshgrid(xlin, ylin)
    cml = generate_cmap(['black', 'black', 'black'])
    cmf = generate_cmap(['darkblue', 'white', 'brown'])
    z = freq[:, i]
    Z = griddata((x, y), z, (X, Y), method='linear')
    plt.ylim(min(y), max(y))
    plt.xlim(min(x), max(x))
    plt.axis('scaled')
    plt.subplot(2, 4, i+1+pn*4)
    interval = np.arange(0, 10, 1)
    plt.title("band #"+str(i), fontsize=9)
    plt.yticks(yt, " ")
    plt.xticks(xt, " ")
    # plt.contour(X, Y, Z, interval, cmap=cml)
    interval = np.arange(0, 11, 1)
    plt.axis('off')
    plt.pcolor(X, Y, Z, vmin=0, vmax=10, cmap=cm.rainbow)
    # plt.colorbar(ticks=[])


def parse_qpoints(qpfile, M):
    f = h5py.File(qpfile)
    freq = f["frequency"]
    qpt = f["qpoint"]
    cqpt = np.dot(qpt, np.transpose(np.linalg.inv((M))))
    if 'A-M' in qpfile:
        y = cqpt[:, 2]
        x = cqpt[:, 1]
    if 'A-K' in qpfile:
        y = cqpt[:, 2]
        x = (cqpt[:, 0]**2 + cqpt[:, 1]**2)**(0.5)*(-1.0)
    return(x, y, freq)


def caserun(qf_A_M, qf_A_K, A, pn):
    x_A_M, y_A_M, freq_A_M = parse_qpoints(qf_A_M, A)
    x_A_K, y_A_K, freq_A_K = parse_qpoints(qf_A_K, A)
    x = np.concatenate((x_A_M, x_A_K), axis=0)
    print x.shape
    print max(x), min(x)
    y = np.concatenate((y_A_M, y_A_K), axis=0)
    print y.shape
    freq = np.concatenate((freq_A_M, freq_A_K), axis=0)

    xt = np.arange(-0.08, 0.08, 0.02)
    yt = np.arange(0.00, 0.18, 0.02)
    for i in range(0, 4):
        contour(i, freq, x, y, xt, yt, pn)
    plt.subplots_adjust(wspace=0.15, hspace=0.15)


def run():
    cA = np.array([[7.8079969443386954, 0.0000000000000000, 0.0000000000000001], [-3.9039984721693477, 6.7619237064685818, 0], [0, 0, 5.6590777347249741]])
    cqpfile_A_M = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/qpoints_A-M.hdf5"
    cqpfile_A_K = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/qpoints_A-K.hdf5"
    sA = np.array([[7.6595137795552795, 0.0000000000000000, 0.0000000000000001], [-3.8297568897776397, 6.6333335137318326, 0], [0, 0, 2.9247116510287272]])
    sqpfile_A_M = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/qpoints_A-M.hdf5"
    sqpfile_A_K = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/qpoints_A-K.hdf5"
    plt.figure(figsize=(16, 8))
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    caserun(cqpfile_A_M, cqpfile_A_K, cA, 0)
    caserun(sqpfile_A_M, sqpfile_A_K, sA, 1)


run()
plt.show()

# plt.savefig("myplot.eps")
# plt.savefig("myplot.png", dpi=800)
# plt.savefig("myplot-colorbar.png", dpi=800)
