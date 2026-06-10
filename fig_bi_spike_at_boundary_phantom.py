#!/usr/bin/env python
import numpy as np
import pickle
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec
from skimage import measure
sys.path.append("/home/kazu/denoise")
# fdir, df are in mlfdev61
fdir = '/home/kazu/restormer_rev2_lim/bi3d/restormer_conv3d/for_single/' +\
       'train/full/211/true_edge/nll/gau2ch/ktrand/rev4/'
df = '/home/kazu/desktop/240424/uNID_data_KO/211/' +\
     'bi3d_scratch_rev4_211_partial_phantom_local_with_gt.pkl'
data_names = ['denoisedx2_5models', 'denoised_5models', 'stride155/expt']


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


def get_images(data_names):
    timage = get_timage(data_names)
    #maskl = get_mask(262)
    #timage[:, :, maskl] = 0.
    refim = timage[2]
    denoisedim = timage[1]
    return refim, denoisedim


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


def pixel_edge_paths(mask: np.ndarray):
    from collections import defaultdict
    """
    右端カラムにある 2値マスク（True/1が内部）の境界を、
    画素格子線上の閉じたポリライン群として返す。
    戻り値: [N_i×2] の配列のリスト。各配列は (x, y) 角点列で先頭=末尾。
    座標系: imshow(..., origin='lower') に重ねられるように -0.5 シフト済み。
    """ 
    m = np.pad(mask.astype(np.uint8), 1, constant_values=0)
    H, W = m.shape

    # 値不一致 → エッジ（水平/垂直）
    hor = m[1:, :] != m[:-1, :]     # (H-1, W)  上下差 ⇒ y と y+1 の間
    ver = m[:, 1:] != m[:, :-1]     # (H, W-1)  左右差 ⇒ x と x+1 の間


    adj = defaultdict(list)
    edges = set()

    def add_edge(p, q):
        adj[p].append(q); adj[q].append(p)
        key = (p, q) if p < q else (q, p)
        edges.add(key)

    # 水平エッジ
    for y in range(H-1):
        for x in np.flatnonzero(hor[y]):
            p = (x,   y+1)
            q = (x+1, y+1)
            add_edge(p, q)

    # 垂直エッジ
    for y in range(H):
        for x in np.flatnonzero(ver[y]):
            p = (x+1, y)
            q = (x+1, y+1)
            add_edge(p, q)

    # エッジ連結 → 閉路ポリライン化
    used = set()
    paths = []
    for e in list(edges):
        if e in used:
            continue
        start = e[0]; cur = e[0]; prev = None
        path = [start]
        while True:
            nxt = []
            for nb in adj[cur]:
                key = (cur, nb) if cur < nb else (nb, cur)
                if key not in used:
                    nxt.append((nb, key))
            if not nxt:
                break
            # 2本のうち前の点を避けて進む
            if prev is None:
                nb, key = nxt[0]
            else:
                cand = [c for c in nxt if c[0] != prev]
                nb, key = cand[0] if cand else nxt[0]
            used.add(key)
            path.append(nb)
            prev, cur = cur, nb
            if nb == start:
                break

        if len(path) >= 2 and path[0] == path[-1]:
            P = np.asarray(path, dtype=float)
            # padding補正(-1) + 角点をimshowのピクセル境界に合わせる(-0.5)
            P -= 1.0
            P -= 0.5
            paths.append(P)

    return paths


def on_frame_xy(p, W, H, atol=1e-9):
    return (np.isclose(p[0], -0.5, atol=atol) or
            np.isclose(p[0], W - 0.5, atol=atol) or
            np.isclose(p[1], -0.5, atol=atol) or
            np.isclose(p[1], H - 0.5, atol=atol))


def inside_xy(p, W, H, atol=1e-9):
    return ((p[0] > -0.5 + atol) and (p[0] < W - 0.5 - atol) and
            (p[1] > -0.5 + atol) and (p[1] < H - 0.5 - atol))


def get_bd(msk):
    bd = pixel_edge_paths(msk)[0]
    nb = bd.shape[0]
    wps = set()
    H, W = msk.shape
    for bidx in range(1, nb-1):
        if on_frame_xy(bd[bidx], W, H):
            if inside_xy(bd[bidx-1], W, H) or inside_xy(bd[bidx+1], W, H):
                wps.add(bidx)
    wp = sorted(wps)
    if on_frame_xy(bd[0], W, H):
        _bd = bd[wp[0]:wp[1]+1]
    else:
        _bd = np.vstack((bd[wp[1]:], bd[0:wp[0]+1]))
    return _bd


def compare_images4_2d():
    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(nrows=3, ncols=4, figure=fig, width_ratios=[1, ]*4,
                           hspace=0.3, wspace=0.03)
    axes = []
    for i in range(3):
        row_ax = []
        for j in range(4):
            a = fig.add_subplot(gs[i, j])
            row_ax.append(a)
        axes.append(row_ax)
    mask = get_mask(smallarea=True)
    # get a list of mask boundary positions, bd: [[y[i], x[i]]].
    bd = measure.find_contours(mask == 1, level=0.5)[0]
    # size of square ROIs
    edge = 16
    # x positions whose y positions == bdpys are picked up from mask boundary
    # positions. They are used to set centers of ROIs.
    bdpys = [20, 61, 53]
    # indicies which select one position from possibly plural positions picked
    # up from each element in bdpys.
    bdidids = [0, 0, 1]
    llimss = []
    ulimss = []
    for bdpy, bdidid in zip(bdpys, bdidids):
        # Note that bd[:, 0] is y value.
        bdids = np.where(bd[:, 0] == bdpy)[0][bdidid]
        bdp = [int(bd[bdids, 1]), int(bd[bdids, 0])]
        llims = [bdp[1]-edge, bdp[0]-edge]
        ulims = [bdp[1]+edge, bdp[0]+edge]
        if llims[0] < 0:
            diff = llims[0]
            llims[0] -= diff
            ulims[0] -= diff
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
    with open(fdir + 'params_scratch_rev4.pkl', 'rb') as f:
        paramimage = pickle.load(f)
    ref = paramimage[:6]
    _, denoised = get_images(data_names)
    for lidx, (ulims, llims) in enumerate(zip(ulimss, llimss)):
        true = ref[4, llims[0]:ulims[0], llims[1]:ulims[1]]
        denoise = denoised[4, llims[0]:ulims[0], llims[1]:ulims[1]]
        _mask = ~mask[llims[0]:ulims[0], llims[1]:ulims[1]]
        vmin = min(np.min(true[true > 0.]), np.min(denoise[denoise > 0.]))
        vmax = max(np.max(true), np.max(denoise))
        ims = []
        for iidx, (data, clabel) in enumerate(
                zip([ref[4], denoised[4], np.abs(ref[4]-denoised[4]), mask],
                    [r'$\mathrm{d_{110}\ /\ \AA}$',
                     r'$\mathrm{d_{110}\ /\ \AA}$',
                     r'$\Delta\mathrm{d_{110}\ /\ \AA}$', 'mask value'])):
            if iidx < 2:
                im, axes[lidx][iidx] = imshow_slice(data, llims, ulims,
                                                    ax=axes[lidx][iidx],
                                                    vmin=vmin, vmax=vmax)
            else:
                im, axes[lidx][iidx] = imshow_slice(data, llims, ulims,
                                                    ax=axes[lidx][iidx])
            ims.append(im)
            map_pos = axes[lidx][iidx].get_position()  # figure 座標
            cax = fig.add_axes([
                map_pos.x1 + 0.002,  # left
                map_pos.y0,                 # bottom（マップと一致）
                0.005,             # width
                map_pos.height              # height（マップと完全一致）
            ])
            cbar = fig.colorbar(im, cax=cax)
            cbar.set_label(clabel, labelpad=0.2)
            cax.tick_params(direction='in', labelsize=8, length=3)
            cbar.locator = MaxNLocator(nbins=4)
            cbar.update_ticks()
        _bd = get_bd(_mask)
        for aidx, _ax in enumerate(axes[lidx][:-1]):
            _ax.plot(_bd[:, 0]+llims[1], _bd[:, 1]+llims[0], ls='-', c='w')
        for ridx in range(3):
            for cidx in range(4):
                axes[ridx][cidx].set_ylabel('y / ch', labelpad=0.1)
                axes[ridx][cidx].set_xlabel('x / ch', labelpad=0.1)
                axes[ridx][cidx].tick_params(length=2, pad=0.01)
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
    #compare_images4_d()
    compare_images4_2d()
