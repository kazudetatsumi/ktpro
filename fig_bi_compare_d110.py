#!/usr/bin/env python
import numpy as np
import os
import pickle
import sys
import copy
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib as mpl
sys.path.append("/home/kazu/denoise")


# Make Arial the default font for all text
mpl.rcParams['font.family'] = 'Arial'         # primary family
mpl.rcParams['font.sans-serif'] = ['Arial']   # optional: ensure sans-serif points to Arial

# Optional: math text style (see section 4)
mpl.rcParams['mathtext.fontset'] = 'dejavusans'  # keeps math readable with a sans style
tfs = '14'
lefs = '10'
tkfs = '12'

mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 9, "axes.labelsize": 11, "xtick.labelsize": 9, "ytick.labelsize": 9,
    "legend.fontsize": 9,
})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'


def get_mask(smallarea=False):
    with open('params_scratch_rev4.pkl', 'rb') as f:
        paramimage = pickle.load(f)
    if not smallarea:
        return paramimage[4] == 0.
    else:
        with open('/home/kazu/desktop/240424/uNID_data_KO/211/bi3d_scratch_rev4_211_partial_phantom_local_with_gt.pkl', 'rb') as f:
            sample_noisy = pickle.load(f)[:, 2:-2, 2:-2]
            openbeam_noisy = pickle.load(f)[:, 2:-2, 2:-2]
            sample_gt = pickle.load(f)[:, 2:-2, 2:-2]
            openbeam_gt = pickle.load(f)[:, 2:-2, 2:-2]
        transmission=sample_gt/openbeam_gt
        sumtransmission=transmission.sum(axis=0)
        bmask = sumtransmission  > 132
        return bmask.T


def get_fitparams(fname='RITS32030_16_40us_Sz1.out'):
    return np.genfromtxt(fname)


def read_flnames(flname='temp_edge_list.dat'):
    print("CHK", flname)
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
    mask = get_mask()
    pos = read_flnames(flname=flname)
    image = np.zeros((len(param_name)*2, mask.shape[0], mask.shape[1]))
    for fidx in range(pos.shape[0]):
        if maskconsider:
            if not mask[pos[fidx, 1], pos[fidx, 0]]:
                image[:, pos[fidx, 1], pos[fidx, 0]] = resultfinal[fidx, :-1]
        else:
            image[:, pos[fidx, 1], pos[fidx, 0]] = resultfinal[fidx, :-1]
    return image


def compare_images4():
    fdir = os.getcwd() + "/"
    data_name = ['denoisedx2_5models', 'denoised_5models', 'stride155/expt']
    for sidx, string in enumerate(data_name):
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
    param_name = ['a0', 'b0', 'a_hkl', 'b_hkl', 'd$_{110}$ / $\\mathrm{\\AA}$',
                  'sigma0']
    mask = get_mask(smallarea=True)
    with open('params_scratch_rev4.pkl', 'rb') as f:
        paramimage = pickle.load(f)
    paramimage[:, mask] = 0.
    timage[:, :, mask] = 0.
    cond2 = np.array(np.where(np.invert(mask)))
    llims = [np.min(cond2[0]), np.min(cond2[1])]
    ulims = [np.max(cond2[0])+1, np.max(cond2[1])+1]
    timage = timage[:, :, llims[0]:ulims[0], llims[1]:ulims[1]]
    mask = mask[llims[0]:ulims[0], llims[1]:ulims[1]]
    x_test = paramimage[:6, llims[0]:ulims[0], llims[1]:ulims[1]]
    x_test_noisy = timage[2]
    bitest = timage[1]
    bitest[6:] = np.array(((timage[1, :6] - timage[0, :6])**2
                          + timage[1, 6:]**2)**0.5)
    cond = ((bitest[4] - bitest[10]) < x_test[4]) & ((bitest[4] + bitest[10])
                                                     > x_test[4])
    print('denoised OK:', np.sum(cond), '/',  np.sum(np.invert(mask)),
          np.sum(cond)/np.sum(np.invert(mask)))
    cond = ((x_test_noisy[4] - x_test_noisy[10]) < x_test[4]) &\
           ((x_test_noisy[4] + x_test_noisy[10]) > x_test[4])
    print('expt OK:', np.sum(cond), '/', np.sum(np.invert(mask)),
          np.sum(cond)/np.sum(np.invert(mask)))
    n = len(param_name)
    r = 4
    plote(x_test_noisy[r], x_test[r], bitest[r], x_test_noisy[r+n],
          np.zeros(x_test[r].shape), bitest[r+n],
          labels=['Observed', 'Ground truth', 'Denoised'])


def plote(
    data, target, pred, datae, targete, prede, savefig='dummy',
    time_idx=[0.2, 0.5, 0.8], q_idx=[0.2, 0.5, 0.8], param_name='Param',
    labels=['Partial', 'Full', 'Partial_denoised'],
):
    fig, ax = plt.subplots(3, 3, figsize=(7.25, 6))
    vmaxim = np.max(
        [np.unique(np.sort(target))[-25], np.unique(np.sort(data))[-25],
         np.unique(np.sort(pred))[-25]])
    vminim = np.min(
            [np.unique(np.sort(target))[15], np.unique(np.sort(data))[15],
             np.unique(np.sort(pred))[15]])
    for didx, (d, label) in enumerate(zip([target, data, pred],
                                          [labels[1], labels[0], labels[2]])):
        _d = copy.deepcopy(d)
        _d[_d == 0.] = np.nan
        cmap = copy.copy(plt.get_cmap('gray_r'))
        cmap.set_bad(color='0.9')
        im = ax[0, didx].imshow(_d, cmap=cmap, origin='lower',
                                interpolation='none', vmin=vminim, vmax=vmaxim)
        ax[0, didx].set_title(label)
        ax[0, didx].set_ylabel('y / ch')
        ax[0, didx].set_xlabel('x / ch')
        ax[0, didx].tick_params(direction='in', axis='both')
    time_idx = [int(t*data.shape[0]) for t in time_idx]
    q_idx = [int(t*data.shape[1]) for t in q_idx]
    tidx = time_idx[1]
    qidx = q_idx[1]
    cond = np.where(data[tidx, :qidx] > 0.)[0]
    ymin1 = np.min(data[tidx, :qidx][cond]-datae[tidx, :qidx][cond])-0.0001
    ymax1 = np.max(data[tidx, :qidx][cond]+datae[tidx, :qidx][cond])+0.0001
    cond = np.where(data[:tidx, qidx] > 0.)[0]
    ymin2 = np.min(data[:tidx, qidx][cond]-datae[:tidx, qidx][cond])-0.0001
    ymax2 = np.max(data[:tidx, qidx][cond]+datae[:tidx, qidx][cond])+0.0001
    margin = 3
    dshape = target.shape
    for (d, label, c, ls) in zip([data[tidx, :qidx], target[tidx, :qidx],
                                  pred[tidx, :qidx]], labels,
                                 ['darkgray', 'k', 'k'], ['-', '--', '-']):
        cond = np.where(d > 0.)[0]
        ax[1, 0].plot(cond, d[d > 0.], label=label, c=c, ls=ls)
        ax[1, 0].set_xlim([cond[0]-margin, cond[-1]+margin])
        ax[1, 0].set_ylim([ymin1, ymax1])
    for didx, (d, de, label, c) in enumerate(
            zip([data[tidx, :qidx], pred[tidx, :qidx]], [datae[tidx, :qidx],
                prede[tidx, :qidx]], [labels[0], labels[2]], ['gray', 'k'])):
        cond = np.where(target[tidx] > 0.)[0]
        ax[1, didx+1].plot(cond, target[tidx][cond], label=labels[1], c='k',
                           ls='--')
        cond = np.where(d > 0.)[0]
        ax[1, didx+1].errorbar(cond, d[cond], yerr=de[cond], marker=".", ms=2,
                               elinewidth=0.3, lw=0, capsize=1.1, label=label,
                               zorder=5, c=c)
        ax[1, didx+1].set_ylim([ymin1, ymax1])
        ax[1, didx+1].set_xlim([cond[0]-margin, cond[-1]+margin])
    ax[0, 0].axhline(tidx, xmin=(cond[0]*1.)/dshape[1],
                     xmax=(cond[-1]*1.)/dshape[1], color='k',
                     linestyle='-', lw=2, alpha=0.8)
    ax[0, 0].axhline(tidx, xmin=(cond[0]*1.)/dshape[1],
                     xmax=(cond[-1]*1.)/dshape[1], color='w',
                     linestyle='-', lw=1)
    margin = 1
    for (d, label, c, ls) in zip([data[:tidx, qidx], target[:tidx, qidx],
                                 pred[:tidx, qidx]],
                                 labels, ['darkgray', 'k', 'k'],
                                 ['-', '--', '-']):
        cond = np.where(d > 0.)[0]
        ax[2, 0].plot(cond, d[cond], label=label, c=c, ls=ls)
        ax[2, 0].set_xlim([cond[0]-margin, cond[-1]+margin])
        ax[2, 0].set_ylim([ymin2, ymax2])
    for didx, (d, de, label, c) in enumerate(
            zip([data[:tidx, qidx], pred[:tidx, qidx]], [datae[:tidx, qidx],
                prede[:tidx, qidx]], [labels[0], labels[2]], ['gray', 'k'])):
        cond = np.where(target[:tidx, qidx] > 0.)[0]
        ax[2, didx+1].plot(cond, target[:tidx, qidx][cond], label=labels[1],
                           c='k', ls='--')
        cond = np.where(d > 0.)[0]
        ax[2, didx+1].errorbar(cond, d[cond], yerr=de[cond], marker=".", ms=2,
                               elinewidth=0.3, lw=0, capsize=1.1, label=label,
                               zorder=5, c=c)
        ax[2, didx+1].set_xlim([cond[0]-margin, cond[-1]+margin])
        ax[2, didx+1].set_ylim([ymin2, ymax2])
    ax[0, 0].axvline(qidx, ymin=(cond[0]*1.)/dshape[0],
                     ymax=(cond[-1]*1.)/dshape[0], color='k',
                     linestyle='-', lw=2, alpha=0.8)
    ax[0, 0].axvline(qidx, ymin=(cond[0]*1.)/dshape[0],
                     ymax=(cond[-1]*1.)/dshape[0], color='w',
                     linestyle='-', lw=1)
    for ridx in range(1, 3):
        for cidx in range(3):
            ax[ridx, cidx].set_ylabel('d$_{110}$ / $\\mathrm{\\AA}$')
            ax[ridx, cidx].tick_params(axis='both', direction='in')
            if ridx == 1:
                ax[ridx, cidx].set_xlabel('x / ch')
            else:
                ax[ridx, cidx].set_xlabel('y / ch')

    plt.tight_layout()
    clabel = r'$\mathrm{d_{110}\ /\ \AA}$'
    map_pos = ax[0, -1].get_position()
    cax = fig.add_axes([
        map_pos.x1 + 0.002,  # left
        map_pos.y0,                 # bottom（マップと一致）
        0.005,             # width
        map_pos.height              # height（マップと完全一致）
    ])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(clabel, labelpad=0.1)
    cax.tick_params(direction='out', labelsize=8, length=2)
    cbar.locator = MaxNLocator(nbins=4)
    cbar.update_ticks()
    plt.show()


compare_images4()
