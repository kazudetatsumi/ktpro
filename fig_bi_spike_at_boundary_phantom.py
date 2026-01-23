#!/usr/bin/env python
import numpy as np
import pickle
import sys
import matplotlib.pyplot as plt
sys.path.append("/home/kazu/denoise")
# fdir, df are in mlfdev61
fdir = '/home/kazu/restormer_rev2_lim/bi3d/restormer_conv3d/for_single/' +\
       'train/full/211/true_edge/nll/gau2ch/ktrand/rev4/'
df = '/home/kazu/desktop/240424/uNID_data_KO/211/' +\
     'bi3d_scratch_rev4_211_partial_phantom_local_with_gt.pkl'


def get_mask(smallarea=False):
    if not smallarea:
        with open(fdir + 'params_scratch_rev4.pkl', 'rb') as f:
            paramimage = pickle.load(f)
        return paramimage[4] == 0.
    else:
        with open(df, 'rb') as f:
            pickle.load(f)
            pickle.load(f)
            sample_gt = pickle.load(f)[:, 2:-2, 2:-2]
            openbeam_gt = pickle.load(f)[:, 2:-2, 2:-2]
        transmission = sample_gt/openbeam_gt
        sumtransmission = transmission.sum(axis=0)
        bmask = sumtransmission > 132
        return bmask.T


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


def imshow_slice(A, llims, ulims, si=1, sj=1, origin='upper',
                 dx=1.0, dy=1.0, x0=0.0, y0=0.0,
                 show_ticks=False, ax=None, **imshow_kwargs):
    """
    A[i0:i1:si, j0:j1:sj] を imshow し、軸を元インデックス（または実座標）で表示。
    実座標: x = j*dx + x0, y = i*dy + y0

    Parameters
    ----------
    A : 2D ndarray
    i0, i1, si : int
    j0, j1, sj : int
    origin : 'upper' or 'lower'
    dx, dy, x0, y0 : float
        実座標のスケーリング（デフォルトはインデックス座標）
    show_ticks : bool
        Trueなら元インデックス/実座標で目盛りを打つ
    ax : matplotlib.axes.Axes or None
    **imshow_kwargs : 追加のimshow引数（cmap等）
    """
    i0 = llims[0]
    i1 = ulims[0]
    j0 = llims[1]
    j1 = ulims[1]
    B = A[i0:i1:si, j0:j1:sj]
    By, Bx = B.shape
    print(B.shape, i0, i1, j0, j1, si, sj)

    # X端
    x_min = (j0 - 0.5*sj) * dx + x0
    x_max = (j0 + sj*(Bx - 1) + 0.5*sj) * dx + x0

    # Y端（originにより上下が入れ替わる）
    if origin == 'upper':
        y_min = (i0 + si*(By - 1) + 0.5*si) * dy + y0
        y_max = (i0 - 0.5*si) * dy + y0
    else:
        y_min = (i0 - 0.5*si) * dy + y0
        y_max = (i0 + si*(By - 1) + 0.5*si) * dy + y0

    if ax is None:
        fig, ax = plt.subplots()

    im = ax.imshow(B, origin=origin,
                   extent=(x_min, x_max, y_min, y_max),
                   interpolation='nearest',
                   **imshow_kwargs)
    if show_ticks:
        xticks = (j0 + np.arange(Bx)*sj) * dx + x0
        yticks = (i0 + np.arange(By)*si) * dy + y0
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
    return im, ax


def compare_images4_2d():
    from skimage import measure
    timage = get_timage(['denoisedx2_5models', 'denoised_5models',
                        'stride155/expt'])
    mask = get_mask(smallarea=True)
    bd = measure.find_contours(mask == 1, level=0.5)[0]
    edge = 16
    bdpys = [20, 61, 53]
    bdidids = [0, 0, 1]
    llimss = []
    ulimss = []
    for bdpy, bdidid in zip(bdpys, bdidids):
        # Note that bd[:, 0] is y value.
        bdids = np.where(bd[:, 0] == bdpy)[0][bdidid]
        bdp = [int(bd[bdids, 1]), int(bd[bdids, 0])]
        llims = [bdp[1]-edge, bdp[0]-edge]
        ulims = [bdp[1]+edge, bdp[0]+edge]
        if llims[1] < 0:
            diff = llims[1]
            llims[1] -= diff
            ulims[1] -= diff
        if ulims[0] > mask.shape[0]:
            diff = ulims[0] - mask.shape[0]
            llims[0] -= diff
            ulims[0] -= diff
        llimss.append(llims)
        ulimss.append(ulims)
    with open('params_scratch_rev4.pkl', 'rb') as f:
        paramimage = pickle.load(f)
    x_test = paramimage[:6]
    bitest = timage[1]
    fig, ax = plt.subplots(3, 4)
    for lidx, (ulims, llims) in enumerate(zip(ulimss, llimss)):
        true = x_test[4, llims[0]:ulims[0], llims[1]:ulims[1]]
        denoise = bitest[4, llims[0]:ulims[0], llims[1]:ulims[1]]
        _mask = ~mask[llims[0]:ulims[0], llims[1]:ulims[1]]
        vmin = min(np.min(true[true > 0.]), np.min(denoise[denoise > 0.]))
        vmax = max(np.max(true), np.max(denoise))
        ax[lidx, 0].imshow(true, vmin=vmin, vmax=vmax)
        ax[lidx, 1].imshow(denoise, vmin=vmin, vmax=vmax)
        ax[lidx, 2].imshow(np.abs(true-denoise))
        ax[lidx, 3].imshow(_mask)
        _bds = measure.find_contours(_mask == 1, level=0.5)[0]
        for aidx, _ax in enumerate(ax[lidx, :-1]):
            _ax.plot(_bds[:, 1], _bds[:, 0], ls=':', c='w')
    plt.show()


def compare_images4_d():
    from scipy import ndimage as ndi
    timage = get_timage(['denoisedx2_5models', 'denoised_5models',
                        'stride155/expt'])
    mask = get_mask(smallarea=True)
    with open('params_scratch_rev4.pkl', 'rb') as f:
        paramimage = pickle.load(f)
    x_test = paramimage[:6]
    bitest = timage[1]
    dist_inside = ndi.distance_transform_edt(~mask)
    dist_outside = ndi.distance_transform_edt(mask)
    abs_errors = np.abs(bitest[4]-x_test[4])
    plt.scatter(dist_inside[dist_inside > 0.].flatten(),
                abs_errors[dist_inside > 0.].flatten(), s=1)
    plt.scatter(dist_outside[dist_outside > 0].flatten()*(-1),
                abs_errors[dist_outside > 0].flatten(), s=1)
    plt.xlim([-10, 25])
    plt.ylabel('Absolute Errors in d$_{110}$ / $\\mathrm{\\AA}$')
    plt.xlabel('Distance from Mask Boundary / ch')
    plt.show()


if __name__ == '__main__':
    compare_images4_d()
    compare_images4_2d()
