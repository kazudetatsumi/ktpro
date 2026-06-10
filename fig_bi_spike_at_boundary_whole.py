#!/usr/bin/env python
import numpy as np
import pickle
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec
from skimage import measure
import cv2
import copy
sys.path.append("/home/kazu/denoise")
# fdir, df are in mlfdev61
fdir = '/home/kazu/restormer_rev2_lim/bi3d/restormer_conv3d/for_single/' +\
       'train/full/211/true_edge/nll/gau2ch/ktrand/'
df = '/home/kazu/desktop/240424/uNID_data_KO/211/openbeam.pkl'
data_names = ['stride_155/expt_nomask', 'denoised_5models_nomask']


with open(df, 'rb') as f:
    openbeam = pickle.load(f)
    openbeam_noisy = pickle.load(f)
    sample = pickle.load(f)
    sample_noisy = pickle.load(f)
    etc = pickle.load(f)
tofini = etc['inilambda']
toffin = etc['finlambda']
dtof = etc['difflambda']
x = np.arange(sample.shape[0])*dtof + tofini
transmission = sample/openbeam*315715/553690
__transmission = transmission.transpose((2, 1, 0))[:, np.newaxis]\
    .astype('float32')


def get_mask(lim):
    return cv2.GaussianBlur(__transmission.squeeze().sum(axis=-1), (5, 5), 0
                            ) > lim*315715/553690.


def get_fitparams(fname='RITS32030_16_40us_Sz1.out'):
    return np.genfromtxt(fname)


def read_flnames(flname='temp_edge_list.dat'):
    pos = []
    for line in open(flname):
        values = line[:-1].split('.')[0].split('_')
        pos.append([int(values[-2]), int(values[-1])])
    return np.array(pos)


def get_paramimage4(fname, fname2, fname3, fname4, flname, maskconsider=False):
    param_name = ['a0', 'b0', 'a_hkl', 'b_hkl', 'd_hkl', 'sigma0']
    result = get_fitparams(fname=fname)
    result2 = get_fitparams(fname=fname2)
    result3 = get_fitparams(fname=fname3)
    result4 = get_fitparams(fname=fname4)
    resultfinal = np.zeros_like(result)
    for icol in range(result.shape[0]):
        resarray = np.array([result[icol, -1], result2[icol, -1],
                             result3[icol, -1], result4[icol, -1]])
        sidx = np.argsort(resarray)
        for nidx, idx in enumerate(sidx):
            if resarray[idx] > 0.001:
                fidx = idx
                break
            elif nidx == len(sidx)-1:
                fidx = idx
        if fidx == 0:
            resultfinal[icol] = result[icol]
        elif fidx == 1:
            resultfinal[icol] = np.insert(np.insert(result2[icol], 5, 3.),
                                          11, 0.)
        elif fidx == 2:
            resultfinal[icol] = np.insert(np.insert(result3[icol], 5, 6.),
                                          11, 0.)
        elif fidx == 3:
            resultfinal[icol] = np.insert(np.insert(result4[icol], 5, 18.),
                                          11, 0.)
    mask = get_mask(262)
    pos = read_flnames(flname=flname)
    image = np.zeros((len(param_name)*2, mask.shape[0], mask.shape[1]))
    for fidx in range(pos.shape[0]):
        if maskconsider:
            if not mask[pos[fidx, 1], pos[fidx, 0]]:
                image[:, pos[fidx, 1], pos[fidx, 0]] = resultfinal[fidx, :-1]
        else:
            image[:, pos[fidx, 1], pos[fidx, 0]] = resultfinal[fidx, :-1]
    return image


def get_timage(data_names):
    for sidx, string in enumerate(data_names):
        fname = fdir + string + '/edge_3.out'
        fname2 = fdir + string + '/edge_3f.out'
        fname3 = fdir + string + '/edge_6f.out'
        fname4 = fdir + string + '/edge_18f.out'
        flname = fdir + string + '/temp_edge_list.dat'
        image = get_paramimage4(fname, fname2, fname3, fname4, flname)
        if sidx == 0:
            timage = image[np.newaxis]
        else:
            timage = np.vstack((timage, image[np.newaxis]))
    return timage


def get_images(data_names):
    timage = get_timage(data_names)
    maskl = get_mask(262)
    timage[:, :, maskl] = 0.
    refim = timage[0]
    denoisedim = timage[1]
    return refim, denoisedim, maskl


def compare_images4_2d():
    fig, ax = plt.subplots(1, 1, figsize=(6, 3))
    mask = get_mask(236)
    bd = measure.find_contours(mask == 1, level=0.5)[0]
    refim, denoisedim, maskl = get_images(data_names)
    data = np.abs(refim[4] - denoisedim[4])
    vmax = np.max(data)
    vmin = np.min(data)
    clabel = r'$\Delta\mathrm{d_{110}\ /\ \AA}$'
    _data = copy.deepcopy(data)
    _data[maskl] = np.nan
    cmap = copy.copy(plt.cm.get_cmap('gray_r'))
    cmap.set_bad(color='0.9')
    im = ax.imshow(_data, cmap=cmap, vmin=vmin, vmax=vmax)
    map_pos = ax.get_position()
    cax = fig.add_axes([
        map_pos.x1 + 0.003,  # left
        map_pos.y0,                 # bottom（マップと一致）
        0.010,             # width
        map_pos.height              # height（マップと完全一致）
    ])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(clabel, labelpad=0.1)
    cax.tick_params(direction='out', labelsize=8, length=2)
    cbar.locator = MaxNLocator(nbins=4)
    cbar.update_ticks()
    ax.plot(bd[:, 1], bd[:, 0], ls='-', c='k', alpha=0.8, lw=2.0)
    ax.plot(bd[:, 1], bd[:, 0], ls='-', c='white', lw=1.0)
    ax.set_ylabel('y / ch', labelpad=0.2)
    ax.set_xlabel('x / ch', labelpad=0.2)
    ax.tick_params(direction='in', length=4)
    plt.show()


if __name__ == '__main__':
    compare_images4_2d()
