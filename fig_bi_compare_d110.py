#!/usr/bin/env python
import numpy as np
import os
import pickle
import sys
import matplotlib.pyplot as plt
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
    param_name = ['a0', 'b0', 'a_hkl', 'b_hkl', 'd$_{110}$ / $\\mathrm{\\AA}$', 'sigma0']
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
    bitest[6:] = np.array(((timage[1, :6] - timage[0, :6])**2 + timage[1, 6:]**2)**0.5)
    cond = ((bitest[4] - bitest[10]) < x_test[4]) & ((bitest[4] + bitest[10]) > x_test[4] )
    print('denoised OK:', np.sum(cond),'/',  np.sum(np.invert(mask)), np.sum(cond)/np.sum(np.invert(mask)))
    cond = ((x_test_noisy[4] - x_test_noisy[10]) < x_test[4]) & ((x_test_noisy[4] + x_test_noisy[10]) > x_test[4] )
    print('expt OK:', np.sum(cond), '/', np.sum(np.invert(mask)), np.sum(cond)/np.sum(np.invert(mask)))
    #showOrigDec(x_test_noisy[:, :, :], x_test[:, :, :], bitest[:, :, :], 
    #            modelname='_restormer_torch_denoise_211_5models_witherror_smallarea_wcbar',
    #            dataname='biparams', num=len(param_name),
    #            param_name=param_name, mask=mask[:, :], labels=['Observed', 'Ground truth', 'Denoised'],
    #            #q_idx=[0.25, 0.55, 0.85], time_idx=[0.2, 0.5, 0.8],
    #            errors=True)
    n = len(param_name)
    r = 4
    plote(x_test_noisy[r], x_test[r], bitest[r], x_test_noisy[r+n],
          np.zeros(x_test[r].shape), bitest[r+n],
          labels=['Observed', 'Ground truth', 'Denoised'])


def showOrigDec(noise, orig, denoise, num=2, pngnum=None, modelname='ae',
                dataname='nr2d', islog=True, param_name=None, mask=None,
                wfile=False, savetxt=False, tint=False,
                labels=['Partial', 'Full', 'Partial_denoised'],
                q_idx=[0.2, 0.5, 0.8], time_idx=[0.2, 0.5, 0.8], errors=False):
    print("ERRORS:", errors)
    if not pngnum:
        pngnum = num
    for i in range(num):
        _noise = noise[i].squeeze()
        _orig = orig[i].squeeze()
        _denoise = denoise[i].squeeze()
        if dataname == 'biparams' and errors:
            if orig.shape[0] == noise.shape[0]/2:
                orig = np.vstack((orig, np.zeros(orig.shape)))
            _noisee = noise[i+len(param_name)].squeeze()
            _orige = orig[i+len(param_name)].squeeze()
            _denoisee = denoise[i+len(param_name)].squeeze()
        if i < pngnum:
            if dataname == 'biparams' and errors is not True:
                print("LINEPlot without errors")
                plot(_noise, _orig, _denoise,
                              savefig='test_' + modelname + '_' + dataname +
                              '-' + str(i)+'.png', param_name=param_name[i],
                              labels=labels, q_idx=q_idx, time_idx=time_idx)
            elif dataname == 'biparams' and errors is True:
                print("LINEPlot with errors")
                plote(_noise, _orig, _denoise, _noisee, _orige, _denoisee,
                               savefig='test_' + modelname + '_' + dataname +
                               '-' + str(i)+'.png', param_name=param_name[i],
                               labels=labels, q_idx=q_idx, time_idx=time_idx)


def plote(data, target, pred, datae, targete, prede, savefig='dummy', time_idx=[0.2, 0.5, 0.8],
         q_idx=[0.2, 0.5, 0.8], param_name='Param',
         labels=['Partial', 'Full', 'Partial_denoised']):
    fig, ax = plt.subplots(3, 3, figsize=(10, 10))
    vmax = np.max([np.max(data[:, 20:]), np.max(target[:, 20:]), np.max(pred[:, 20:])])
    vmin = np.min([np.unique(np.sort(data))[1], np.unique(np.sort(target))[1],
                   np.unique(np.sort(pred))[1]])
    vmax = np.max([np.unique(np.sort(target))[-10], np.unique(np.sort(pred))[-30]])
    vmin = np.min([np.unique(np.sort(target))[10],
                   np.unique(np.sort(pred))[15]])
    vmax = np.max([np.unique(np.sort(target))[-15], np.unique(np.sort(data+datae))[-15],
                   np.unique(np.sort(pred))[-15]])
    vmin = np.min([np.unique(np.sort(target))[5], np.unique(np.sort(data-datae))[5],
                   np.unique(np.sort(pred))[5]])
    vmaxim = np.max([np.unique(np.sort(target))[-25], np.unique(np.sort(data))[-25],
                   np.unique(np.sort(pred))[-25]])
    vminim = np.min([np.unique(np.sort(target))[15], np.unique(np.sort(data))[15],
                   np.unique(np.sort(pred))[15]])
    #vmax = 2.045
    #vmin = 2.02
    print(f'vmax={vmax}')
    print(f'vmin={vmin}')
    ax[0, 0].imshow(data, interpolation='none', vmin=vminim, vmax=vmaxim)
    #ax[0, 0].imshow(data, interpolation='none')
    ax[0, 0].set_title(labels[0], fontsize=tfs)
    ax[0, 0].set_ylabel('y / ch', fontsize=tkfs)
    ax[0, 0].set_xlabel('x / ch', fontsize=tkfs)
    ax[0, 0].tick_params(axis='both', labelsize=lefs)
    im = ax[0, 2].imshow(target,  interpolation='none', vmin=vminim, vmax=vmaxim)
    #ax[0, 1].imshow(target,  interpolation='none')
    ax[0, 2].set_yticks([])
    ax[0, 2].set_title(labels[1], fontsize=tfs)
    ax[0, 2].set_xlabel('x / ch', fontsize=tkfs)
    ax[0, 2].tick_params(axis='both', labelsize=lefs)
    ax[0, 1].imshow(pred, interpolation='none', vmin=vminim, vmax=vmaxim)
    #ax[0, 2].imshow(pred, interpolation='none')
    ax[0, 1].set_yticks([])
    ax[0, 1].set_title(labels[2], fontsize=tfs)
    ax[0, 1].set_xlabel('x / ch', fontsize=tkfs)
    ax[0, 1].tick_params(axis='both', labelsize=lefs)
    #cbar = fig.colorbar(im, ax=ax[0, 2], fraction=0.046, pad=0.04)
    #cbar.ax.tick_params(labelsize=12)
    #cbar.set_label(param_name, fontsize=tfs)


    N = 3
    time_idx = [int(t*data.shape[0]) for t in time_idx]
    q_idx = [int(t*data.shape[1]) for t in q_idx]
    for i in range(0, N):
        #ax[1, i].plot(data[time_idx[i]], 'o', ms=3, label=labels[0])
        #ax[1, i].errorbar(range(data[time_idx[i]].shape[0]), data[time_idx[i]],
        #                  yerr=datae[time_idx[i]], marker="o", ms=3,
        #                  elinewidth=0.4, lw=0, capsize=2, label=labels[0], zorder=0)
        #ax[1, i].plot(pred[time_idx[i]], label=labels[2])
        ax[1, i].errorbar(range(pred[time_idx[i]].shape[0]), pred[time_idx[i]],
                          yerr=prede[time_idx[i]], marker="o", ms=3,
                          elinewidth=0.4, lw=0, capsize=3, label=labels[2], zorder=5)
        #ax[1, i].errorbar(range(pred[time_idx[i]].shape[0]), pred[time_idx[i]],
        #                  yerr=np.abs((data[time_idx[i]]-pred[time_idx[i]])**2+prede[time_idx[i]]**2)**0.5, marker="o", ms=3,
        #                  elinewidth=0.4, lw=0, capsize=3, label=labels[2], zorder=5)
        ax[1, i].plot(target[time_idx[i]], label=labels[1], zorder=10)
        #ax[1, i].errorbar(range(target[time_idx[i]].shape[0]), target[time_idx[i]],
        #                  yerr=targete[time_idx[i]], marker="o", ms=3,
        #                  elinewidth=1, lw=0, capsize=3, label=labels[1])
        ax[1, i].set_xlabel('x / ch', fontsize=tfs)
        leg = ax[1, i].legend(frameon=False, fontsize=tkfs, handlelength=1.2,
                        title=f'at y = {time_idx[i]}', alignment='left')
        leg.set_title(f'at y= {time_idx[i]}', prop={'size': 12})
        ax[1, i].set_ylim([vmin, vmax])
        ax[1, i].tick_params(axis='both', labelsize=tkfs)
        ax[0, 0].axhline(time_idx[i], color='r', linestyle='--', lw=1)
        if i > 0:
            ax[1, i].set_yticks([])
        #ax[2, i].plot(np.array(data[:, q_idx[i]]), 'o', ms=3, label=labels[0])
        #ax[2, i].errorbar(range(data[:, q_idx[i]].shape[0]), data[:, q_idx[i]],
        #                  yerr=datae[:, q_idx[i]], marker="o", ms=2,
        #                  elinewidth=0.4, lw=0, capsize=2, label=labels[0])
        #ax[2, i].plot(np.array(pred[:, q_idx[i]]), label=labels[2])
        ax[2, i].errorbar(range(pred[:, q_idx[i]].shape[0]), pred[:, q_idx[i]],
                          yerr=prede[:, q_idx[i]], marker="o", ms=3,
                          elinewidth=0.4, lw=0, capsize=3, label=labels[2])
        #ax[2, i].errorbar(range(pred[:, q_idx[i]].shape[0]), pred[:, q_idx[i]],
        #                  yerr=np.abs((data[:, q_idx[i]] - pred[:, q_idx[i]])**2+prede[:, q_idx[i]]**2)**0.5, marker="o", ms=3,
        #                  elinewidth=0.4, lw=0, capsize=3, label=labels[2])
        ax[2, i].plot(np.array(target[:, q_idx[i]]), label=labels[1])
        #ax[2, i].errorbar(range(target[:, q_idx[i]].shape[0]), target[:, q_idx[i]],
        #                  yerr=targete[:, q_idx[i]], marker="o", ms=3,
        #                  elinewidth=1, lw=0, capsize=3, label=labels[1])
        ax[2, i].set_xlabel('y / ch', fontsize=tfs)
        leg = ax[2, i].legend(frameon=False, fontsize=tkfs, handlelength=1.2,
                        title=f'at x= {q_idx[i]}', alignment='left')
        leg.set_title(f'at x= {q_idx[i]}', prop={'size': 12})
        ax[2, i].set_ylim([vmin, vmax])
        ax[2, i].tick_params(axis='both', labelsize=tkfs)
        ax[0, 0].axvline(q_idx[i], color='r', linestyle='--', lw=1)
        if i > 0:
            ax[2, i].set_yticks([])
    ax[1, 0].set_ylabel(param_name, fontsize=tfs)
    ax[2, 0].set_ylabel(param_name, fontsize=tfs)
    for j in range(1, 3):
        yrange = np.array([ax[j, i].get_ylim() for i in range(0, N)])
        for i in range(0, N):
            ax[j, i].set_ylim(np.min(yrange[:, 0]), np.max(yrange[:, 1]))
    plt.tight_layout()
    if savefig != '':
        if not ('.png' in savefig or '.pdf' in savefig):
            savefig += '.pdf'
        plt.savefig(savefig)
    else:
        plt.show()
    plt.close()


def plot(data, target, pred, savefig='', time_idx=[0.2, 0.5, 0.8],
         q_idx=[0.2, 0.5, 0.8], param_name='Param',
         labels=['Partial', 'Full', 'Partial_denoised']):
    fig, ax = plt.subplots(3, 3, figsize=(10, 10))
    vmax = np.max([np.max(data[:, 20:]), np.max(target[:, 20:]), np.max(pred[:, 20:])])
    vmin = np.min([np.unique(np.sort(data))[1], np.unique(np.sort(target))[1],
                   np.unique(np.sort(pred))[1]])
    vmax = np.max([np.unique(np.sort(target))[-10], np.unique(np.sort(pred))[-10]])
    vmax = np.max([np.unique(np.sort(target))[-10], np.unique(np.sort(pred))[-10], np.unique(np.sort(data))[-130]])
    vmax = np.max([np.unique(np.sort(target))[-10], np.unique(np.sort(pred))[-10]])
    vmin = np.min([np.unique(np.sort(target))[20],
                   np.unique(np.sort(pred))[20]])
    print(f'vmax={vmax}')
    print(f'vmin={vmin}')
    ax[0, 0].imshow(data, interpolation='none', vmin=vmin, vmax=vmax)
    ax[0, 0].set_title(labels[0], fontsize=tfs)
    ax[0, 0].set_ylabel('y / ch', fontsize=tfs)
    ax[0, 0].set_xlabel('x / ch', fontsize=tfs)
    ax[0, 1].imshow(target,  interpolation='none', vmin=vmin, vmax=vmax)
    ax[0, 1].set_yticks([])
    ax[0, 1].set_title(labels[1], fontsize=tfs)
    ax[0, 1].set_xlabel('x / ch', fontsize=tfs)
    im = ax[0, 2].imshow(pred, interpolation='none', vmin=vmin, vmax=vmax)
    ax[0, 2].set_yticks([])
    ax[0, 2].set_title(labels[2], fontsize=tfs)
    ax[0, 2].set_xlabel('x / ch', fontsize=tfs)
    #cbar = fig.colorbar(im, ax=ax[0, 2], fraction=0.046, pad=0.04)
    #cbar.ax.tick_params(labelsize=12)
    #cbar.set_label(param_name, fontsize=tfs)
    N = 3
    time_idx = [int(t*data.shape[0]) for t in time_idx]
    q_idx = [int(t*data.shape[1]) for t in q_idx]
    for i in range(0, N):
        #ax[1, i].plot(data[time_idx[i]], 'o', ms=3, label=labels[0])
        ax[1, i].plot(data[time_idx[i]], label=labels[0])
        ax[1, i].plot(target[time_idx[i]], label=labels[1])
        ax[1, i].plot(pred[time_idx[i]], label=labels[2])
        ax[1, i].set_xlabel('x / ch', fontsize=tfs)
        leg = ax[1, i].legend(frameon=False, fontsize=tkfs, handlelength=1.2,
                        title=f'at y = {time_idx[i]}', alignment='left')
        leg.set_title(f'at y= {time_idx[i]}', prop={'size': 12})
        ax[1, i].set_ylim([vmin, vmax])
        ax[0, 0].axhline(time_idx[i], color='r', linestyle='--', lw=1)
        if i > 0:
            ax[1, i].set_yticks([])
        #ax[2, i].plot(np.array(data[:, q_idx[i]]), 'o', ms=3, label=labels[0])
        ax[2, i].plot(np.array(data[:, q_idx[i]]), ms=3, label=labels[0])
        ax[2, i].plot(np.array(target[:, q_idx[i]]), label=labels[1])
        ax[2, i].plot(np.array(pred[:, q_idx[i]]), label=labels[2])
        ax[2, i].set_xlabel('y / ch', fontsize=tfs)
        leg = ax[2, i].legend(frameon=False, fontsize=tkfs, handlelength=1.2,
                        title=f'at x= {q_idx[i]}', alignment='left')
        leg.set_title(f'at x= {q_idx[i]}', prop={'size': 12})
        ax[2, i].set_ylim([vmin, vmax])
        ax[0, 0].axvline(q_idx[i], color='r', linestyle='--', lw=1)
        if i > 0:
            ax[2, i].set_yticks([])
    ax[1, 0].set_ylabel(param_name, fontsize=tfs)
    ax[2, 0].set_ylabel(param_name, fontsize=tfs)
    for j in range(1, 3):
        yrange = np.array([ax[j, i].get_ylim() for i in range(0, N)])
        for i in range(0, N):
            ax[j, i].set_ylim(np.min(yrange[:, 0]), np.max(yrange[:, 1]))
    plt.tight_layout()
    if savefig != '':
        if not ('.png' in savefig or '.pdf' in savefig):
            savefig += '.pdf'
        plt.savefig(savefig)
    else:
        plt.show()
    plt.close()





compare_images4()
