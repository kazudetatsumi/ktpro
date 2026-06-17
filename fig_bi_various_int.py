#!/usr/bin/env python
import numpy as np
import os
import pickle
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter
sys.path.append("/home/kazu/denoise")
from tools import showOrigDec


def get_fitparams(fname='RITS32030_16_40us_Sz1.out'):
    return np.genfromtxt(fname)


def get_tofIn():
    return np.arange(152)*20+23000


def get_gtparamimage():
    with open('../../../rev4/params_scratch_rev4.pkl', 'rb') as f:
        paramimage = pickle.load(f)
    return paramimage


def read_flnames(flname='temp_edge_list.dat'):
    pos = []
    for line in open(flname):
        values = line[:-1].split('.')[0].split('_')
        pos.append([int(values[-2]), int(values[-1])])
    return np.array(pos)


def get_data(datapath):
    with open(datapath, 'rb') as f:
        sample = pickle.load(f)[:, 2:-2, 2:-2].transpose((2, 1, 0))
        openbeam = pickle.load(f)[:, 2:-2, 2:-2].transpose((2, 1, 0))
    return sample, openbeam


def _get_mask(smallarea=False):
    paramimage = get_gtparamimage()
    if not smallarea:
        return paramimage[4] == 0.
    else:
        points = np.where(paramimage[4] != 0.)
        points = np.array(list(zip(points[0], points[1])))
        center = np.mean(points, axis=0)
        centered_points = points - center
        scaled_points = centered_points * 0.95
        new_points = np.int32(scaled_points + center)
        rows = [coord[0] for coord in new_points]
        cols = [coord[1] for coord in new_points]
        mask = np.full(paramimage[4].shape, True, dtype=bool)
        mask[(rows, cols)] = False
        return mask


def get_mask(smallarea=False):
    paramimage = get_gtparamimage()
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


def evaluate_PSNR(_noise, _orig):
    from skimage import metrics
    val_min = _orig.min()
    val_range = _orig.max() - val_min
    __orig = (_orig - val_min)/val_range
    __noise = (_noise - val_min)/val_range
    return metrics.peak_signal_noise_ratio(__orig, __noise, data_range=1.0)


def get_PSNR(noise, orig):
    psnrs = np.zeros((noise.shape[0]))
    for idx in range(noise.shape[0]):
        psnrs[idx] = evaluate_PSNR(noise[idx], orig[idx])
    return psnrs


def compare_images4(data_name):
    fdir = "/home/kazu/restormer_rev2_lim/bi3d/restormer_conv3d/" +\
           "for_single/train/full/211/true_edge/nll/gau2ch/ktrand/" + \
           "wosingle/stride155/various_int/rev4/"
    fname = fdir + data_name + '/edge_3.out'
    fname2 = fdir + data_name + '/edge_3f.out'
    fname3 = fdir + data_name + '/edge_6f.out'
    fname4 = fdir + data_name + '/edge_18f.out'
    flname = fdir + data_name + '/temp_edge_list.dat'
    image = get_paramimage4(fname, fname2, fname3, fname4, flname)
    param_name = ['a0', 'b0', 'a_hkl', 'b_hkl', 'd_hkl', 'sigma0']
    psnr = np.zeros((1, 2, len(param_name))) # extraxis for vstack, smallarea, params
    paramimage = get_gtparamimage()
    for cidx, (smallarea, tag) in enumerate(zip([False, True], ["", "_smallarea"])):
        mask = get_mask(smallarea=smallarea)
        paramimage[:, mask] = 0.
        image[:, mask] = 0.
        cond2 = np.array(np.where(np.invert(mask)))
        llims = [np.min(cond2[0]), np.min(cond2[1])]
        ulims = [np.max(cond2[0])+1, np.max(cond2[1])+1]
        _image = image[:, llims[0]:ulims[0], llims[1]:ulims[1]]
        _mask = mask[llims[0]:ulims[0], llims[1]:ulims[1]]
        x_test = paramimage[:6, llims[0]:ulims[0], llims[1]:ulims[1]]
        if cidx==1:
            fig, ax = plt.subplots(1, 2)
            ax[0].imshow(_image[4], vmin=2.02, vmax=2.04)
            ax[0].set_title(data_name)
            ax[1].imshow(x_test[4], vmin=2.02, vmax=2.04)
            ax[1].set_title('gt')
            plt.savefig(data_name+'.png')
            plt.close()
        psnr[0, cidx] = get_PSNR(_image[:6], x_test)
    return psnr


def err_transmission(sample, openbeam, fac):
    err_sample = sample**0.5
    err_openbeam = openbeam**0.5
    return ((1./openbeam*err_sample)**2 +
            (sample/openbeam**2*err_openbeam)**2)**0.5*fac


def write_rits_inputfiles_with_error(sample, openbeam, fac, amp):
    mask = get_mask(smallarea=False)
    outDir = "./exptx"+str(amp)[:-2]+'stride155/'
    transmission = sample/openbeam*fac
    x = get_tofIn()
    if not os.path.isdir(outDir):
        os.makedirs(outDir)
    for vpos in range(sample.shape[0]):
        for hpos in range(sample.shape[1]):
            if not mask[vpos, hpos]:
                output = np.vstack((x, transmission[vpos, hpos]))
                output = np.vstack((output,
                                    err_transmission(sample[vpos, hpos],
                                                     openbeam[vpos, hpos],
                                                     fac)
                                    ))
                np.savetxt(outDir+'exptx'+str(amp)[:-2]+'stride155_'+str(hpos)
                           + '_'+str(vpos)+'.txt', output.T)


def samplerun1():
    fac = 315715./553690.
    for amp in 2**np.linspace(0, 10, 11):
        dataprepath = '/home/kazu/desktop/240424/uNID_data_KO/255/' +\
                      'bi3d_scratch_rev4_211_phantom_local_intx'
        datapostpath = '_stride155.pkl'
        datapath = dataprepath + str(amp)[:-2] + datapostpath
        sample, openbeam = get_data(datapath)
        write_rits_inputfiles_with_error(sample, openbeam, fac, amp)


def plotpsnr(psnr):
    lfs = '12'
    fig, ax = plt.subplots(3, 2, figsize=(6, 10))
    sidx = 1
    for pidx, paramname in enumerate(['a0', 'b0', 'a_hkl', 'b_hkl',
                                      'd$_{110}$', 'sigma0']):
        icol = pidx % 2
        irow = pidx // 2
        ax[irow, icol].plot(2**np.linspace(0, 10, 11), psnr[:-11, sidx, pidx],
                            marker='x', c='k', lw=1, ls='--')
        ax[irow, icol].plot([1., 2., 4., 8., 16., 32., 64., 128., 256., 512.,
                             1024.],
                            psnr[-11:, sidx, pidx], marker='o', c='k', lw=1,
                            ls='--')
        ax[irow, icol].set_xscale('log', base=2)
        ax[irow, icol].set_xlabel('Intenisty coefficient, $\\alpha$',
                                  fontsize=lfs)
        ax[irow, icol].set_ylabel('PSNR / dB', fontsize=lfs)
        ax[irow, icol].set_title(paramname)
        ax[irow, icol].tick_params(direction='in', axis='both', which='both', labelsize=lfs)
        # Major ticks: 2^-3, 2^0, 2^3, 2^6, 2^9
        major_ticks = [2**i for i in [0, 3, 6, 9]]
        ax[irow, icol].set_xticks(major_ticks)
        # Minor ticks: 2倍刻み（base=2, subs=[1] → 2^n のみ）
        # subs=[1] だと 2^n のみ、subs=[1,2,3,...] で細かく制御可能
        ax[irow, icol].xaxis.set_minor_locator(LogLocator(base=2, subs=[1.0],
                                               numticks=100))
        ax[irow, icol].xaxis.set_minor_formatter(NullFormatter())
        # Major tickラベルを 2^n 形式に
        labels = [f"$2^{{{i}}}$" for i in [0, 3, 6, 9]]
        ax[irow, icol].set_xticklabels(labels)
    plt.tight_layout()
    #plt.savefig('various_int_phantom_bi3d_psnr_.png')
    plt.show()


def samplerun2():
    for aidx, amp in enumerate(2**np.linspace(0, 10, 11)):
        exptdir = 'exptx'+str(amp)[:-2]+'stride155'
        if aidx == 0:
            _psnr = compare_images4(exptdir)
            psnr = _psnr
        else:
            psnr = np.vstack((psnr, compare_images4(exptdir)))
    for denoiseddir in ['denoised_full', 'denoised_x2', 'denoised_x4',
                        'denoised_x8', 'denoised_x16', 'denoised_x32',
                        'denoised_x64', 'denoised_x128', 'denoised_x256',
                        'denoised_x512', 'denoised_x1024']:
        psnr = np.vstack((psnr, compare_images4(denoiseddir)))
    plotpsnr(psnr)


samplerun2()
