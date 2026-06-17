#!/usr/bin/env python
# Because a possiblly useful code for exlaining models uses pytorch,
# I wrote a denoising code for pytorch.
# Kazuyoshi TATSUMI 2024/04/18
# CUDA_VISIBLE_DEVICES=0,1 torchrun --nnodes 1 --nproc_per_node 2 sample.py
import numpy as np
from tqdm import tqdm
import torch.nn as nn
from torch.utils.data import DataLoader
import torch
import torch.optim as optim
import os
import pickle
import sys
import matplotlib as mpl
from restormer_arch_var8 import Restormer
sys.path.append("/home/kazu/denoise")
import matplotlib.pyplot as plt
#from torch.utils.data.distributed import DistributedSampler
from tools import showOrigDec
#from torchmodels import UNet
from torchtools import EarlyStopping, noisedDataset, train_model_ddp
from explain import perturb_explanation
from tifffile import imread
np.random.seed(314)
torch.manual_seed(314)
torch.cuda.manual_seed(314)
torch.backends.cudnn.deterministic = True
#rank = int(os.environ["LOCAL_RANK"])
#torch.cuda.set_device(rank)
#world_size = torch.cuda.device_count()
#torch.distributed.init_process_group(backend='nccl', init_method='env://',
#                                     world_size=world_size)

mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 9, "axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 9,
})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

def get_tofIn():
    return np.arange(152)*20+23000


def calc_mean(tdata):
    return np.mean(tdata, axis=0)


def calc_var(tvar, tdata):
    M = tdata.shape[0]
    #std = (1/M*(tvar.sum(axis=0) + (tdata**2).sum(axis=0)) - np.mean(tdata, axis=0)**2)**0.5
    return 1./M*(tvar.sum(axis=0) + (tdata**2).sum(axis=0))\
        - np.mean(tdata, axis=0)**2

gathered_inferred_data = 'outputt.pkl'
with open('/home/kazu/desktop/240424/uNID_data_KO/211/bi3d_scratch_rev4_211_partial_phantom_local_with_gt.pkl',
          'rb') as f:
    sample_noisy = pickle.load(f)[:, 2:-2, 2:-2]
    openbeam_noisy = pickle.load(f)[:, 2:-2, 2:-2]
    sample_gt = pickle.load(f)[:, 2:-2, 2:-2]
    openbeam_gt = pickle.load(f)[:, 2:-2, 2:-2]


x = get_tofIn()
transmission_noisy = sample_noisy/openbeam_noisy
__transmission_noisy = transmission_noisy.\
    transpose((2, 1, 0))[:, np.newaxis].astype('float32')
__sample_noisy = sample_noisy.transpose((2, 1, 0))[np.newaxis, np.newaxis]\
    .astype('float32')
__openbeam_noisy = openbeam_noisy.transpose((2, 1, 0))[np.newaxis, np.newaxis]\
    .astype('float32')
__sample_gt = sample_gt.transpose((2, 1, 0))[np.newaxis, np.newaxis]\
    .astype('float32')
__openbeam_gt = openbeam_gt.transpose((2, 1, 0))[np.newaxis, np.newaxis]\
    .astype('float32')
div = 1.0
tsfms = torch.from_numpy
exptset = noisedDataset(__sample_noisy, __openbeam_noisy, tsfms)
batch_size = 1
exptloader = DataLoader(exptset, batch_size=batch_size, num_workers=2,
                        pin_memory=True, shuffle=False)
paths = ['../seed0/chk_54_-0.12946_-0.07817.pt',
         '../seed1/chk_8_-0.09366_-0.05927.pt',
         '../seed2/chk_9_-0.16705_-0.03814.pt',
         '../seed3/chk_9_-0.07207_-0.08073.pt',
         '../seed4/chk_5_-0.09606_-0.12873.pt'
         ]

if not os.path.exists(gathered_inferred_data):
    device = "cuda:0"
    model = Restormer(dim=12, out_channels=2).to(device)
    dirty, clean = next(iter(exptloader))
    dirty, clean = dirty.to(device), clean.to(device)
    print('constructing ' + gathered_inferred_data)
    for pidx, path in enumerate(paths):
        checkpoint = torch.load(path, map_location=torch.device('cpu'))
        model.load_state_dict(checkpoint['model'])
        model.eval()
        with torch.no_grad():
            nrtest, varnrtest = model(dirty)
            onrtest, varonrtest = model(clean)
        nrtest = nrtest.to('cpu').detach().numpy().squeeze()[
                np.newaxis, :, np.newaxis]
        onrtest = onrtest.to('cpu').detach().numpy().squeeze()[
                np.newaxis, :, np.newaxis]
        varnrtest = varnrtest.to('cpu').detach().numpy().squeeze()[
                np.newaxis, :, np.newaxis]
        varonrtest = varonrtest.to('cpu').detach().numpy().squeeze()[
                np.newaxis, :, np.newaxis]
        if pidx == 0:
            nrtestt = nrtest
            onrtestt = onrtest
            varnrtestt = varnrtest
            varonrtestt = varonrtest
        else:
            nrtestt = np.vstack((nrtestt, nrtest))
            onrtestt = np.vstack((onrtestt, onrtest))
            varnrtestt = np.vstack((varnrtestt, varnrtest))
            varonrtestt = np.vstack((varonrtestt, varonrtest))
    varnrtest = calc_var(varnrtestt, nrtestt)
    varonrtest = calc_var(varonrtestt, onrtestt)
    nrtest = calc_mean(nrtestt)
    onrtest = calc_mean(onrtestt)
    with open(gathered_inferred_data, 'wb') as f:
        pickle.dump(nrtest, f, 4)
        pickle.dump(varnrtest, f, 4)
        pickle.dump(onrtest, f, 4)
        pickle.dump(varonrtest, f, 4)
else:
    with open(gathered_inferred_data, 'rb') as f:
        nrtest = pickle.load(f)
        varnrtest = pickle.load(f)
        onrtest = pickle.load(f)
        varonrtest = pickle.load(f)

__sample_noisy = __sample_noisy.squeeze()[:, np.newaxis]
__openbeam_noisy = __openbeam_noisy.squeeze()[:, np.newaxis]
__sample_gt = __sample_gt.squeeze()[:, np.newaxis]
__openbeam_gt = __openbeam_gt.squeeze()[:, np.newaxis]


def plot_profiles():
    vpos = 20 - 16
    hpos = 35 - 8
    fig, ax = plt.subplots(7, 2, figsize=(8, 12))
    vmax = 1.0
    vmin = 0.
    for didx, (denoised, expt, dstr, estr) in enumerate(
            zip([fullnrtest/fullonrtest*315715./553690.],
                [__transmission],
                ['full_denoised'],
                ['full_asr'])):
        ax[0, didx].imshow(expt[vpos, 0, :, :], vmax=vmax, vmin=vmin, aspect=3)
        ax[0, didx].set_title(estr + ' at y=' + str(vpos) + ' ch')
        ax[0, didx].set_xlabel(r'$\lambda$ / ch')
        ax[0, didx].set_ylabel('x / ch')
        ax[1, didx].imshow(denoised[vpos, 0, :, :], vmin=vmin, vmax=vmax,
                           aspect=3)
        ax[1, didx].set_title(dstr + ' at y=' + str(vpos) + ' ch')
        ax[1, didx].set_xlabel(r'$\lambda$ / ch')
        ax[1, didx].set_ylabel('x / ch')
        ax[2, didx].imshow(expt[:, 0, hpos, :]*1.03, vmax=vmax, vmin=vmin,
                           aspect=3)
        ax[2, didx].set_title(estr + ' at x=' + str(hpos) + ' ch')
        ax[2, didx].set_xlabel(r'$\lambda$ / ch')
        ax[2, didx].set_ylabel('y / ch')
        ax[3, didx].imshow(denoised[:, 0, hpos, :], vmax=vmax, vmin=vmin,
                           aspect=3)
        ax[3, didx].set_title(dstr + ' at x=' + str(hpos) + ' ch')
        ax[3, didx].set_xlabel(r'$\lambda$ / ch')
        ax[3, didx].set_ylabel('y / ch')
        ax[4, didx].plot(x, expt[vpos, 0, hpos], label=estr)
        ax[4, didx].plot(x, denoised[vpos, 0, hpos], label=dstr)
        ax[4, didx].set_title(estr+', '+dstr+' at (x,y)=('+str(hpos)+',' +
                              str(vpos)+')'+' ch')
        ax[4, didx].set_xlabel('TOF microsec')
        ax[4, didx].set_ylabel('Transmission')
        ax[4, didx].set_ylim([0.4, 1.1])
        ax[4, didx].legend()
        ax[5, didx].imshow(denoised[:, 0, :, 50], vmin=vmin, vmax=vmax)
        ax[5, didx].set_title(dstr + ' at TOF=50 ch')
        ax[5, didx].set_ylabel('y / ch')
        ax[5, didx].set_xlabel('x / ch')
        ax[6, didx].imshow(expt[:, 0, :, 50], vmin=vmin, vmax=vmax)
        ax[6, didx].set_title(estr + ' at TOF=50 ch')
        ax[6, didx].set_ylabel('y / ch')
        ax[6, didx].set_xlabel('x / ch')
    plt.tight_layout()
    plt.savefig('bi2d_restormer_results.png')
    plt.show()


def get_mask():
    print(__transmission_noisy.shape)
    #sum_transmission = np.sum(__transmission_noisy.squeeze(), axis=2)
    #plt.imshow(sum_transmission > 148)
    #plt.show()
    #print(type(sum_transmission))
    #plt.hist(sum_transmission.flatten(), bins='auto')
    #plt.show()
    #plt.plot(__transmission[5, 0, 20, :])
    #plt.show()
    #___transmission = np.average((nrtest/onrtest)[:, 0, :, 90:], axis=2)
    #___transmission = np.average(__transmission[:, 0, :, 90:], axis=2)
    #plt.imshow(__transmission[:, 0, :, :].sum(axis=2))
    #plt.show()
    #return ___transmission >= 0.95
    #return paramimage[4] == 0.
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


def plot_masked():
    ___transmission = np.average(__transmission[:, 0, :, 345:355], axis=2)
    #mask = ___transmission >= 0.8
    mask = get_mask()
    print("CHCK", mask.shape)
    print(mask[27, 44])
    ___transmission[mask] = 0.
    plt.imshow(___transmission)
    plt.show()


def get_fitparams(fname='RITS32030_16_40us_Sz1.out'):
    #fname = 'fit_results_expt_7240-28000.txt'
    #fname = 'RITS32030_16_40us_Sz1-3.out.expt_noisy'
    #fname = 'RITS32030_16_40us_Sz1-4.out'
    return np.genfromtxt(fname)


def plot_fitparams():
    param_name = ['PRODEN / ' + r'${\rm 10^{22}\ cm^{-2}}$',
                  r'LATCON / ${\rm \AA}$', 'MICRST / um',
                  'MDCOFF', 'CRSIZE / um']
    result = get_fitparams()
    mask = get_mask()
    fig, ax = plt.subplots(5, 1, figsize=(10, 10))
    for pidx, param in enumerate(param_name):
        image = np.zeros((mask.shape[0], mask.shape[1]))
        fidx = 0
        for hpos in range(__sample.shape[2]):
            for vpos in range(__sample.shape[0]):
                if not mask[vpos, hpos]:
                    image[vpos, hpos] = result[fidx, pidx]
                    fidx += 1
        ax[pidx].imshow(image, vmin=np.min(result[:, pidx]),
                        vmax=np.max(result[:, pidx]))
        ax[pidx].set_title(param)
    plt.tight_layout()
    plt.show()


def read_flnames(flname='temp_edge_list.dat'):
    print("CHK", flname)
    pos = []
    for line in open(flname):
        values = line[:-1].split('.')[0].split('_')
        pos.append([int(values[-2]), int(values[-1])])
    return np.array(pos)


def plot_fitparams2():
    param_name = ['PRODEN / ' + r'${\rm 10^{22}\ cm^{-2}}$',
                  r'LATCON / ${\rm \AA}$', 'MICRST / um', 'MDCOFF',
                  'CRSIZE / um']
    result = get_fitparams()
    mask = get_mask()
    pos = read_flnames()
    fig, ax = plt.subplots(5, 1, figsize=(10, 10))
    for pidx, param in enumerate(param_name):
        image = np.zeros((mask.shape[0], mask.shape[1]))
        for fidx in range(pos.shape[0]):
            image[pos[fidx, 1], pos[fidx, 0]] = result[fidx, pidx]
        ax[pidx].imshow(image, vmin=np.min(result[:, pidx]),
                        vmax=np.max(result[:, pidx]))
        ax[pidx].set_title(param)
    plt.tight_layout()
    plt.show()


def get_paramimage(fname, flname):
    param_name = ['PRODEN / ' + r'${\rm 10^{22}\ cm^{-2}}$',
                  r'LATCON / ${\rm \AA}$', 'MICRST / um', 'MDCOFF',
                  'CRSIZE / um']
    result = get_fitparams(fname=fname)
    mask = get_mask()
    pos = read_flnames(flname=flname)
    image = np.zeros((len(param_name), mask.shape[0], mask.shape[1]))
    for pidx, param in enumerate(param_name):
        for fidx in range(pos.shape[0]):
            image[pidx, pos[fidx, 1], pos[fidx, 0]] = result[fidx, pidx]
    return image


def get_paramimage2(fname, fname2, flname):
    param_name = ['PRODEN / ' + r'${\rm 10^{22}\ cm^{-2}}$',
                  r'LATCON / ${\rm \AA}$', 'MICRST / um', 'MDCOFF',
                  'CRSIZE / um']
    #param_name = ['PRODEN / ' + r'${\rm 10^{22}\ cm^{-2}}$',
    #              r'LATCON / ${\rm \AA}$', 'MICRST / um', 'MDCOFF']
    result = get_fitparams(fname=fname)
    result2 = get_fitparams(fname=fname2)
    cond = result[:, -1] > result2[:, -1]
    # if the squared kai error in the case of the 5th parameter fixed as zero
    # is smaller, this result is selected for the 1st, 2nd, 3rd, and 4th
    # parameters and the 5th parameter is set as zero.
    result[:, :4][cond] = result2[:, :4][cond]
    result[:, 4][cond] = 0.
    mask = get_mask()
    pos = read_flnames(flname=flname)
    image = np.zeros((len(param_name), mask.shape[0], mask.shape[1]))
    for pidx, param in enumerate(param_name):
        for fidx in range(pos.shape[0]):
            image[pidx, pos[fidx, 1], pos[fidx, 0]] = result[fidx, pidx]
            #image[pidx, pos[fidx, 1], pos[fidx, 0]] = result2[fidx, pidx]
    return image


def get_paramimage4(fname, fname2, fname3, fname4, flname, maskconsider=False):
    param_name = ['a0', 'b0', 'a_hkl', 'b_hkl', 'd_hkl', 'sigma0']
    #param_name = ['PRODEN / ' + r'${\rm 10^{22}\ cm^{-2}}$',
    #              r'LATCON / ${\rm \AA}$', 'MICRST / um', 'MDCOFF',
    #              'CRSIZE / um']
    print(fname)
    result = get_fitparams(fname=fname)
    result2 = get_fitparams(fname=fname2)
    result3 = get_fitparams(fname=fname3)
    result4 = get_fitparams(fname=fname4)
    print(result.shape) # 275 11
    print(result2.shape) # 275 11
    print(result3.shape) # 275 11
    print(result4.shape) # 275 11
    resultfinal = np.zeros_like(result)
    #for icol in range(result.shape[0]):
    #    if result[icol, -1] < result2[icol, -1] and result[icol, -1] < result3[icol, -1] and result[icol, -1] < result4[icol, -1]:
    #        resultfinal[icol] = result[icol]
    #    if result2[icol, -1] < result[icol, -1] and result2[icol, -1] < result3[icol, -1] and result2[icol, -1] < result3[icol, -1]:
    #        resultfinal[icol] = np.insert(np.insert(result2[icol], 5, 3.), 11, 0.)
    #    if result3[icol, -1] < result2[icol, -1] and result3[icol, -1] < result[icol, -1] and result3[icol, -1] < result4[icol, -1]:
    #        resultfinal[icol] = np.insert(np.insert(result3[icol], 5, 6.), 11, 0.)
    #    if result4[icol, -1] < result2[icol, -1] and result4[icol, -1] < result[icol, -1] and result4[icol, -1] < result3[icol, -1]:
    #        resultfinal[icol] = np.insert(np.insert(result4[icol], 5, 18.), 11, 0.)
    # The above procedure contained a bug that not converged parameters in result is necessarily selected.
    # The following procedure fixes it.
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


def compare_images():
    data_name = ['expt', 'denoise', 'denoise_noisy', 'expt_noisy']
    for sidx, string in enumerate(data_name):
        fname = 'RITS32030_16_40us_Sz1.out.' + string
        flname = 'temp_edge_list.dat.' + string
        image = get_paramimage(fname, flname)
        if sidx == 0:
            timage = image[np.newaxis]
        else:
            timage = np.vstack((timage, image[np.newaxis]))
    param_name = ['PRODEN / ' + r'${\rm 10^{22}\ cm^{-2}}$',
                  r'LATCON / ${\rm \AA}$', 'MICRST / um', 'MDCOFF',
                  'CRSIZE / um']
    fig, ax = plt.subplots(len(param_name), len(data_name), figsize=(10, 10))
    for sidx, string in enumerate(data_name):
        for pidx, param in enumerate(param_name):
            ax[pidx, sidx].imshow(np.abs(timage[sidx, pidx]),
                                  vmin=np.unique(np.sort(np.abs(
                                                 timage[:, pidx])))[1],
                                  vmax=np.max(timage[:, pidx]))
            if pidx == 0:
                ax[pidx, sidx].set_title(string + '\n' + param)
            else:
                ax[pidx, sidx].set_title(param)
            print(np.unique(np.sort(np.abs(timage[:, pidx]))))
    plt.tight_layout()
    plt.show()


def compare_images2():
    fdir = "/home/kazu/denoise/"
    fdir = "/home/kazu/RITS_fit/"
    data_name = ['expt_lim', 'denoise_rev2_lim',
            #'denoise_noisy_rev2_lim_131', 'expt_noisy_lim_divopt']
                 'denoise_noisy_rev2_lim_divopt', 'expt_noisy_lim_divopt']
    for sidx, string in enumerate(data_name):
        fname = fdir +  'RITS32030_16_40us_Sz1.out.' + string
        fname2 = fdir + 'RITS32030_16_40us_Sz0.out.' + string
        flname = fdir + 'temp_edge_list.dat.' + string
        image = get_paramimage2(fname, fname2, flname)
        if sidx == 0:
            timage = image[np.newaxis]
        else:
            timage = np.vstack((timage, image[np.newaxis]))
    param_name = ['PRODEN / ' + r'${\rm 10^{22}\ cm^{-2}}$',
                  r'LATCON / ${\rm \AA}$', 'MICRST / um', 'MDCOFF',
                  'CRSIZE / um']
    #param_name = ['PRODEN / ' + r'${\rm 10^{22}\ cm^{-2}}$',
    #              r'LATCON / ${\rm \AA}$', 'MICRST / um', 'MDCOFF']
    mask = get_mask()
    cond2 = np.array(np.where(np.invert(mask)))
    llims = [np.min(cond2[0]), np.min(cond2[1])]
    ulims = [np.max(cond2[0])+1, np.max(cond2[1])+1]
    ## select the smallest rectangle area containing the unmasked area.
    timage = timage[:, :, llims[0]:ulims[0], llims[1]:ulims[1]]
    mask = mask[llims[0]:ulims[0], llims[1]:ulims[1]]
    x_test = np.abs(timage[0])
    x_test_noisy = np.abs(timage[3])
    bitest = np.abs(timage[2])
    #bitest = np.abs(timage[1])
    showOrigDec(x_test, x_test_noisy, bitest,
                modelname='restormer_torch_denoise_rev3_divopt', dataname='biparams',
                num=len(param_name), islog=False, param_name=param_name,
                mask=mask)
    fig, ax = plt.subplots(len(param_name), len(data_name), figsize=(10, 6))
    for sidx, string in enumerate(data_name):
        for pidx, param in enumerate(param_name):
            ax[pidx, sidx].imshow(np.abs(timage[sidx, pidx]),
                                  vmin=np.unique(np.sort(np.abs(
                                                 timage[:, pidx])))[1],
                                  vmax=np.max(timage[:, pidx]))
            if pidx == 0:
                ax[pidx, sidx].set_title(string + '\n' + param)
            else:
                ax[pidx, sidx].set_title(param)
    plt.tight_layout()
    plt.show()

def compare_images4():
    fdir = os.getcwd() + "/" 
    #fdir = "/home/kazu/RITS_fit/"
    data_name = ['expt_noisy', 'denoised_noisy_5models']
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
    param_name = ['a0', 'b0', 'a_hkl', 'b_hkl', 'd_hkl', 'sigma0']
    mask = get_mask()
    with open('params_scratch_rev.pkl', 'rb') as f:
        paramimage = pickle.load(f)
    print(paramimage.shape)
    print(timage[0].shape)
    paramimage[:, mask] = 0.
    timage[:, :, mask] = 0.
    #plt.imshow(mask)
    #plt.show()
    cond2 = np.array(np.where(np.invert(mask)))
    llims = [np.min(cond2[0]), np.min(cond2[1])]
    ulims = [np.max(cond2[0])+1, np.max(cond2[1])+1]
    timage = timage[:, :, llims[0]:ulims[0], llims[1]:ulims[1]]
    mask = mask[llims[0]:ulims[0], llims[1]:ulims[1]]
    x_test = paramimage[:6, llims[0]:ulims[0], llims[1]:ulims[1]]
    x_test_noisy = timage[0]
    #x_test = timage[1]
    bitest = timage[1]
    print('chk', x_test.shape)
    print('chk', x_test_noisy.shape)
    print('chk', bitest.shape)
    print('chk', mask.shape)
    showOrigDec(x_test_noisy[:, :, :], x_test[:, :, :], bitest[:, :, :], 
                modelname='restormer_torch_denoise_211_5models_smallarea',
                dataname='biparams', num=len(param_name), islog=False,
                param_name=param_name, mask=mask[:, :], labels=['expt_noisy', 'gt', 'denoised_noisy'],
                errors=False)
    fig, ax = plt.subplots(len(param_name), len(data_name), figsize=(30, 30))
    for sidx, string in enumerate(data_name):
        for pidx, param in enumerate(param_name):
            ax[pidx, sidx].imshow(np.abs(timage[sidx, pidx]),
                                  vmin=np.unique(np.sort(np.abs(
                                      timage[:, pidx, :, 20:])))[1],
                                  vmax=np.max(timage[:, pidx, :, 20:]))
            if pidx == 0:
                ax[pidx, sidx].set_title(string.split("_lim_strd")[0] +
                                         '\n' + param)
            else:
                ax[pidx, sidx].set_title(param)
    plt.tight_layout()
    plt.savefig('bi_params.png')
    plt.show()


def plot_energyfiltered():
    vmax = 1.0
    vmin = 0.
    fig, ax = plt.subplots(2, 2, figsize=(10, 10))
    for didx, (denoised, expt, dstr, estr) in enumerate(
            zip([nrtest/onrtest,
                 fullnrtest/fullonrtest*315715./553690.],
                [__transmission_noisy, __transmission],
                ['partial_denoised', 'full_denoised'],
                ['partial_asr', 'full_asr'])):
        ax[0, didx].imshow(denoised[:, 0, :, 350], vmin=vmin, vmax=vmax)
        ax[0, didx].set_title(dstr + ' at TOF=350 ch')
        ax[0, didx].set_ylabel('y / ch')
        ax[0, didx].set_xlabel('x / ch')
        ax[1, didx].imshow(expt[:, 0, :, 350], vmin=vmin, vmax=vmax)
        ax[1, didx].set_title(estr + ' at TOF=350 ch')
        ax[1, didx].set_ylabel('y / ch')
        ax[1, didx].set_xlabel('x / ch')
    plt.tight_layout()
    plt.show()


def plot_bi2d():
    vmax = 1.0
    vmin = 0.
    vpos = 20 - 16
    hpos = 35 - 8
    #vpos = 12
    fig, ax = plt.subplots(4, 2, figsize=(10, 10))
    for didx, (denoised, expt, dstr, estr) in enumerate(
            zip([nrtest/onrtest,
                 fullnrtest/fullonrtest*315715./553690.],
                [__transmission_noisy, __transmission],
                ['partial_denoised', 'full_denoised'],
                ['partial_asr', 'full_asr'])):
        ax[0, didx].imshow(denoised[vpos, 0, :, :], vmin=vmin, vmax=vmax, aspect=3)
        ax[0, didx].set_title('transmittance ' + dstr + ' at y=' + str(vpos) + ' ch')
        ax[0, didx].set_ylabel('x / ch')
        ax[0, didx].set_xlabel('t / ch')
        ax[1, didx].imshow(expt[vpos, 0, :, :], vmin=vmin, vmax=vmax, aspect=3)
        ax[1, didx].set_title('transmittance ' + estr + ' at y=' + str(vpos) + ' ch')
        ax[1, didx].set_ylabel('x / ch')
        ax[1, didx].set_xlabel('t / ch')
        ax[2, didx].imshow(denoised[:, 0, hpos, :], vmin=vmin, vmax=vmax, aspect=6)
        ax[2, didx].set_title('transmittance ' + dstr + ' at x=' + str(hpos) + ' ch')
        ax[2, didx].set_ylabel('y / ch')
        ax[2, didx].set_xlabel('t / ch')
        ax[3, didx].imshow(expt[:, 0, hpos, :], vmin=vmin, vmax=vmax, aspect=6)
        ax[3, didx].set_title('transmittance ' + estr + ' at x=' + str(hpos) + ' ch')
        ax[3, didx].set_ylabel('y / ch')
        ax[3, didx].set_xlabel('t / ch')
    plt.tight_layout()
    plt.savefig(f'bi2dspectra_denoise_{hpos}_{vpos}.png')
    plt.show()


def plot_spectra():
    vpos = 30
    hpos = 100
    #vpos = 12
    fig, ax = plt.subplots(4, 1, figsize=(10, 10))
    for didx, (denoised, expt, dstr, estr) in enumerate(
            zip([nrtest,
                 fullnrtest,
                 onrtest,
                 fullonrtest],
                 [__sample_noisy, __sample/div,
                 __openbeam_noisy, __openbeam/div],
                ['partial_denoised', 'full_denoised',
                 'partial_denoised_openbeam', 'full_denoised_openbeam'],
                ['partial_asr', 'full_asr',
                 'partial_openbeam_asr', 'full_openbeam_asr'])):
        print(expt.shape)
        ax[didx].set_title(f'BI {dstr} and {estr} at ({hpos}, {vpos}))')
        ax[didx].scatter(x, expt[vpos, 0, hpos], marker='o', s=1.5)
        ax[didx].plot(x, denoised[vpos, 0, hpos])
    plt.tight_layout()
    plt.savefig(f'bi_spectra_denoise_{hpos}_{vpos}.png')
    plt.show()
    ylim = [np.max((__sample_noisy/__openbeam_noisy)[vpos, 0, hpos])-0.45, np.max((__sample_noisy/__openbeam_noisy)[vpos, 0, hpos])]
    fig, ax = plt.subplots(2, 1, figsize=(10, 10))
    for didx, (denoised, expt, dstr, estr) in enumerate(
            zip([nrtest/onrtest,
                 fullnrtest/fullonrtest*315715./553690],
                [__sample_noisy/__openbeam_noisy,
                 __sample/__openbeam*315715./553690],
                ['partial_denoised', 'full_denoised'],
                ['partial_asr', 'full_asr'])
            ):
        print(denoised.shape)
        ax[didx].scatter(x, expt[vpos, 0, hpos], marker='o', s=1.5)
        ax[didx].plot(x, denoised[vpos, 0, hpos])
        ax[didx].set_title(f'transmission {dstr} and {estr}  at ({hpos}, {vpos})')
        ax[didx].set_ylabel('transmission')
        ax[didx].set_xlabel('t / ch')
        #ax[didx].set_ylim(ylim)
    plt.tight_layout()
    plt.savefig(f'transmission_spectra_denoise_{hpos}_{vpos}.png')
    plt.show()


def plot_spectra_with_error():
    vpos = 30
    hpos = 100 
    #vpos = 12
    fig, ax = plt.subplots(1, 2, figsize=(10, 10))
    for didx, (denoised, expt, gt, dstr, estr, gstr) in enumerate(
            zip([np.vstack((nrtest.transpose((1, 0, 2, 3)), varnrtest.transpose((1, 0, 2, 3)))),
                 np.vstack((onrtest.transpose((1, 0, 2, 3)), varonrtest.transpose((1, 0, 2, 3))))],
                [__sample_noisy, __openbeam_noisy],
                [__sample_gt, __openbeam_gt],
                ['full_denoised', 'full_denoised_openbeam'],
                ['full_asr', 'full_openbeam_asr'],
                ['full_asr_gt', 'full_openbeam_asr_gt'])):
        print(expt.shape, denoised.shape, gt.shape)
        ax[didx].set_title(f'BI {dstr} and {estr} at ({hpos}, {vpos}))')
        ax[didx].scatter(x, expt[vpos, 0, hpos], marker='o', s=1.5)
        ax[didx].plot(x, denoised[0, vpos, hpos])
        ax[didx].plot(x, gt[vpos, 0, hpos], c='r')
        ax[didx].fill_between(x, denoised[0, vpos, hpos] - denoised[1, vpos, hpos]**0.5 , denoised[0, vpos, hpos] + denoised[1, vpos, hpos]**0.5, alpha=0.3, label='Uncertainty')
    plt.tight_layout()
    plt.savefig(f'bi_spectra_denoise_{hpos}_{vpos}.png')
    plt.show()
    divvar = (nrtest/onrtest)**2 * (varnrtest/(nrtest**2) + varonrtest/(onrtest**2))
    fig, ax = plt.subplots(1, 1, figsize=(10, 20))
    for didx, (denoised, expt, gt, dstr, estr, gstr) in enumerate(
            zip([np.vstack(((nrtest/onrtest).transpose((1, 0, 2, 3)), divvar.transpose((1, 0, 2, 3))))],
                [__sample_noisy/__openbeam_noisy],
                [__sample_gt/__openbeam_gt],
                ['full_denoised'],
                ['full_asr'],
                ['full_asr_gt'])
            ):
        print(denoised.shape)
        ax.scatter(x, expt[vpos, 0, hpos], marker='o', s=1.5)
        ax.plot(x, denoised[0, vpos, hpos])
        ax.plot(x, gt[vpos, 0, hpos], c='r')
        ax.fill_between(x, denoised[0, vpos, hpos] - denoised[1, vpos, hpos]**0.5 , denoised[0, vpos, hpos] + denoised[1, vpos, hpos]**0.5, alpha=0.3, label='Uncertainty')
        ax.set_title(f'Transmission {dstr} and {estr}  at ({hpos}, {vpos})', fontsize=18)
        ax.set_ylabel('Transmission', fontsize=18)
        ax.set_xlabel('t / ch', fontsize=18)
    plt.tight_layout()
    plt.savefig(f'transmission_spectra_denoise_{hpos}_{vpos}.png')
    plt.show()


def plot_spectra_with_error_united():
    vpos = 30
    hpos = 100
    divvar = (nrtest/onrtest)**2 * (varnrtest/(nrtest**2) +
                                    varonrtest/(onrtest**2))
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    for didx, (denoised, expt, gt, dstr, estr, gstr, vstr) in enumerate(
            zip([np.vstack((nrtest.transpose((1, 0, 2, 3)),
                            varnrtest.transpose((1, 0, 2, 3)))),
                 np.vstack((onrtest.transpose((1, 0, 2, 3)),
                            varonrtest.transpose((1, 0, 2, 3)))),
                 np.vstack(((nrtest/onrtest).transpose((1, 0, 2, 3)),
                            divvar.transpose((1, 0, 2, 3))))],
                [__sample_noisy, __openbeam_noisy,
                    __sample_noisy/__openbeam_noisy],
                [__sample_gt, __openbeam_gt, __sample_gt/__openbeam_gt],
                ['denoised', 'denoised_openbeam', 'denoised'],
                ['$I$', '$I_0$', '$I/I_{0}$'],
                ['asr_gt', 'openbeam_asr_gt', 'asr_gt'],
                ['Neutron count', 'Neutron count', 'Transmission'])):
        if didx == 0:
            print(expt.shape, denoised.shape, gt.shape)
            #ax.set_title(f'{estr}')
            ax.scatter(x, expt[vpos, 0, hpos], marker='o', s=5, c='k')
            ax.fill_between(
                x, denoised[0, vpos, hpos] - denoised[1, vpos, hpos]**0.5,
                denoised[0, vpos, hpos] + denoised[1, vpos, hpos]**0.5, alpha=1.0,
                color='gray', label='Uncertainty')
            ax.plot(x, denoised[0, vpos, hpos], c='w', lw=1.3)
            ax.plot(x, gt[vpos, 0, hpos], c='k', ls='--', lw=1.3)
            ax.set_xlabel(r'tof / $\mu$s')
            ax.set_ylabel(vstr)
            ax.tick_params(axis='both', direction='in')
    plt.tight_layout()
    plt.savefig('fig_bi_spectra_denoise.eps')
    plt.show()


def err_transmission(sample, openbeam, fac):
    err_sample = sample**0.5
    err_openbeam = openbeam**0.5
    return ((1./openbeam*err_sample)**2 +
            (sample/openbeam**2*err_openbeam)**2)**0.5*fac


def res(coeff):
    mask = get_mask()
    t = __transmission.squeeze()[np.invert(mask)]
    t_noisy = __transmission_noisy.squeeze()[np.invert(mask)]
    return (t - t_noisy*coeff).flatten()


def optimize(variable=0.999):
    import scipy.optimize as so
    out = so.least_squares(res, variable)
    print(out.x[0])
    return out.x[0]


def write_rits_inputfiles():
    mask = get_mask()
    plt.imshow(mask)
    plt.show()
    #const = optimize()
    #chk_6_38.01_40.32.ptt
    outDir = "./2D-BI-OUT_6_3801_4032/"
    if not os.path.isdir(outDir):
        os.makedirs(outDir)
    for vpos in range(__sample.shape[0]):
        for hpos in range(__sample.shape[2]):
            if not mask[vpos, hpos]:
                #output_noisy = np.vstack((x,
                #                          __transmission_noisy[vpos, 0, hpos,
                #                                               :]*const))
                #output_noisy = np.vstack((output_noisy,
                #                          err_transmission(__sample_noisy[
                #                                           vpos, 0, hpos, :],
                #                                           __openbeam_noisy[
                #                                           vpos, 0, hpos, :],
                #                                           1.)
                #                          ))
                #np.savetxt(outDir+'expt_noisy_'+str(hpos)+'_'+str(vpos)+'.txt',
                #           output_noisy.T)
                output = np.vstack((x, __transmission[vpos, 0, hpos]))
                output = np.vstack((output,
                                    err_transmission(__sample[vpos, 0, hpos],
                                                     __openbeam[vpos, 0, hpos],
                                                     315715./553690)
                                    ))
                np.savetxt(outDir+'expt_phantom_'+str(hpos)+'_'+str(vpos)+'.txt',
                           output.T)
                #output_noisy_denoised = np.vstack((x,
                #                                  (nrtest /
                #                                      onrtest)
                #                                   [vpos, 0, hpos]*const))
                #output_noisy_denoised = np.vstack((output_noisy_denoised,
                #                                  err_transmission(
                #                                      __sample_noisy[vpos, 0,
                #                                                     hpos],
                #                                      __openbeam_noisy[vpos, 0,
                #                                                       hpos],
                #                                      1.)
                #                                   ))
                #output_noisy_denoised = np.vstack((output_noisy_denoised,
                #                                   np.zeros((x[7240//40:28000//40+1].shape)
                #                                            ) + 0.01
                #                                   ))
                #np.savetxt(outDir+'denoised_noisy_experr_'+str(hpos)+'_' +
                #           str(vpos) + '.txt', output_noisy_denoised.T)
                output_denoised = np.vstack((x,
                                            (fullnrtest /
                                             fullonrtest)[vpos, 0,
                                             hpos] * 315715/553690))
                output_denoised = np.vstack((output_denoised,
                                             err_transmission(
                                                 __sample[vpos, 0, hpos],
                                                 __openbeam[vpos, 0, hpos],
                                                 315715./553690)
                                             ))
                #output_denoised = np.vstack((output_denoised,
                #                             np.zeros((x[7240//40:28000//40+1].shape))
                #                             + 0.01))
                np.savetxt(outDir + 'denoised_phantom_experr_' + str(hpos) + '_' +
                           str(vpos) + '.txt', output_denoised.T)


def write_rits_inputfiles_with_error():
    mask = get_mask()
    plt.imshow(mask)
    plt.show()
    outDir = "./2D-BI-OUT_5models/"
    divvar = varnrtest/onrtest**2 + varonrtest*(nrtest**2)/(onrtest**4)
    if not os.path.isdir(outDir):
        os.makedirs(outDir)
    for vpos in range(__sample_noisy.shape[0]):
        for hpos in range(__sample_noisy.shape[2]):
            if not mask[vpos, hpos]:
                output = np.vstack((x, __transmission_noisy[vpos, 0, hpos]))
                output = np.vstack((output,
                                    err_transmission(__sample_noisy[vpos, 0, hpos],
                                                     __openbeam_noisy[vpos, 0, hpos],
                                                     1.)
                                    ))
                np.savetxt(outDir+'expt_noisy_'+str(hpos)+'_'+str(vpos)+'.txt',
                           output.T)
                #output = np.vstack((x, __transmission[vpos, 0, hpos]))
                #output = np.vstack((output,
                #                    err_transmission(__sample[vpos, 0, hpos],
                #                                     __openbeam[vpos, 0, hpos],
                #                                     315715./553690)
                #                    ))
                #np.savetxt(outDir+'expt_'+str(hpos)+'_'+str(vpos)+'.txt',
                #           output.T)
                output_denoised = np.vstack((x, (nrtest/onrtest)[vpos, 0, hpos]))
                output_denoised = np.vstack((output_denoised,
                                             divvar[vpos, 0, hpos]**0.5
                                             ))
                np.savetxt(outDir + 'denoised_noisy_' + str(hpos) + '_' +
                           str(vpos) + '.txt', output_denoised.T)


def reuse_images4_for_phantom_data(string):
    fdir = os.getcwd() + "/" + string + "/"
    fname = fdir + 'edge_3.out'
    fname2 = fdir + 'edge_3f.out'
    fname3 = fdir + 'edge_6f.out'
    fname4 = fdir + 'edge_18f.out'
    flname = fdir + 'temp_edge_list.dat'
    with open('paramimage_5models.pkl', 'wb') as f:
        pickle.dump(get_paramimage4(fname, fname2, fname3, fname4, flname),
                    f, 4)


#plot_masked()
#plot_fitparams()
#compare_images()
#compare_images2()
#compare_images4()
#write_rits_inputfiles()
#write_rits_inputfiles_with_error()
#optimize()
#plot_bi2d()
#plot_energyfiltered()
#plot_spectra()
#plot_spectra_with_error()
plot_spectra_with_error_united()
#plot_profiles()
#reuse_images4_for_phantom_data('expt_5models')
