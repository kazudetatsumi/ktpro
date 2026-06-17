#!?usr/bin/env python
import pickle
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.gridspec import GridSpec
import numpy as np
import os
from scipy.ndimage import binary_dilation, distance_transform_edt, gaussian_filter1d
#plt.rcParams["font.family"] = "Arial"   # 使用するフォント
plt.rcParams["font.size"] = 10
plt.rcParams['mathtext.default'] = 'regular'



from scipy.ndimage import gaussian_filter, sobel
np.random.seed(314)


def get_mask(expt=True):
    if expt:
        sample_expt, openbeam_expt, transmission_expt = get_data()
        import cv2
        return cv2.GaussianBlur(transmission_expt.sum(axis=-1), (5, 5), 0) > 236*315715./553690
    else:
        sample_phantom, openbeam_phantom, transmission_phantom_gt = get_data(expt=False)
        sumtransmission_phantom = transmission_phantom_gt.sum(axis=-1)
        bmask = sumtransmission_phantom > 132*315715./553690
        return bmask


def get_data(expt=True):
    if expt:
        with open('/home/kazu/desktop/240424/uNID_data_KO/211/openbeam.pkl',
                  'rb') as f:
            openbeam_expt = pickle.load(f)
            openbeam_expt_noisy = pickle.load(f)
            sample_expt = pickle.load(f)
            sample_expt_noisy = pickle.load(f)
            etc = pickle.load(f)
        sample_expt = sample_expt.transpose((2, 1, 0))
        openbeam_expt = openbeam_expt.transpose((2, 1, 0))
        transmission = sample_expt/openbeam_expt*315715./553690
        return sample_expt, openbeam_expt, transmission
    else:
        #with open('/home/kazu/desktop/240424/uNID_data_KO/211/bi3d_scratch_rev4_211_partial_phantom_local_with_gt.pkl', 'rb') as f:
        with open('/home/kazu/desktop/240424/uNID_data_KO/211/bi3d_scratch_rev4_211_phantom_local.pkl', 'rb') as f:
            sample_noisy = pickle.load(f)[:, 2:-2, 2:-2]
            openbeam_noisy = pickle.load(f)[:, 2:-2, 2:-2]
            #sample_gt = pickle.load(f)[:, 2:-2, 2:-2]
            #openbeam_gt = pickle.load(f)[:, 2:-2, 2:-2]
        with open('/home/kazu/desktop/240424/uNID_data_KO/211/bi3d_scratch_rev4_211_partial_phantom_local_with_gt.pkl', 'rb') as f:
            _sample_noisy = pickle.load(f)[:, 2:-2, 2:-2]
            _openbeam_noisy = pickle.load(f)[:, 2:-2, 2:-2]
            sample_gt = pickle.load(f)[:, 2:-2, 2:-2]
            openbeam_gt = pickle.load(f)[:, 2:-2, 2:-2]
        sample_gt = sample_gt.transpose((2, 1, 0))
        sample_noisy = sample_noisy.transpose((2, 1, 0))
        openbeam_noisy = openbeam_noisy.transpose((2, 1, 0))
        openbeam_gt = openbeam_gt.transpose((2, 1, 0))
        transmission_gt = sample_gt/openbeam_gt*315715./553690
        return sample_noisy, openbeam_noisy, transmission_gt


def get_simudata():
    gathered_data = "/home/kazu/3D-BI/bi3d_testbcc_simudata_rev2_lim_single_resize_full_211_true_edgev_ktrand.pkl"
    with open(gathered_data, 'rb') as f:
        sample = pickle.load(f)
    num_data = sample.shape[0]
    #noisy = np.random.poisson(sample)
    train_size = int(num_data*0.8)
    train_idx = np.random.choice(num_data, train_size,  replace=False)
    x_train = sample[train_idx]
    print("random.poisson is being  appyied now")
    x_train_noisy = np.random.poisson(x_train[:])
    return x_train_noisy.squeeze()


def grad_map_xy(T3d, mask=None, sigma=0.5, agg='max'):
    # T3d: (H, W, T)  ; mask3d optional bool
    #plt.plot(T3d[30, 100])
    #plt.show()
    if sigma > 0:
        T3d = gaussian_filter(T3d, sigma=(sigma, sigma, 0))
    gx = sobel(T3d, axis=0, mode='reflect')
    gy = sobel(T3d, axis=1, mode='reflect')
    g = np.hypot(gx, gy)  # (H, W, T)
    if agg == 'max':
        g2d = g.max(axis=2)
    elif agg == 'p95':
        g2d = np.quantile(g, 0.95, axis=2)
    else:
        raise ValueError
    #g2d = g[:, :, 120:].sum(axis=-1)
    #g2d[mask] = 0.
    #import copy
    #_g2d = copy.deepcopy(g2d)
    #if mask is not None:
    #    _g2d[mask] = 0.
    #plt.imshow(_g2d)
    #plt.show()
    #if mask is not None:
    #    #m2d = mask3d.any(axis=2)
    #    return g2d[np.invert(mask)]
    if mask is not None:
        return g2d[np.invert(mask)].ravel()
    else:
        return g2d.ravel()


def grad_map_2d(T3d, mask=None, sigma=0.5, agg='max'):
    # T3d: (H, W, T)  ; mask3d optional bool
    #plt.plot(T3d[30, 100])
    #plt.show()
    if sigma > 0:
        T3d = gaussian_filter(T3d, sigma=(sigma, sigma, 0))
    gx = sobel(T3d, axis=0, mode='reflect')
    gy = sobel(T3d, axis=1, mode='reflect')
    g = np.hypot(gx, gy)  # (H, W, T)
    if agg == 'max':
        g2d = g.max(axis=2)
    elif agg == 'p95':
        g2d = np.quantile(g, 0.95, axis=2)
    else:
        raise ValueError
    if mask is not None:
        g2d[mask] = 0.
    return g2d


def hist_density(values, bins):
    hist, edges = np.histogram(values, bins=bins, density=True)
    centers = 0.5*(edges[:-1] + edges[1:])
    return centers, hist


def survival_curve(g, grid=None):
    g = np.sort(g[np.isfinite(g)])
    if grid is None:
        grid = np.quantile(g, np.linspace(0, 0.999, 200))
    # S(g)=P(G>g) を階段状に評価
    S = 1.0 - np.searchsorted(g, grid, side='right')/g.size
    return grid, S


def compare_grads():
    gsfile = 'gs.pkl'
    if not os.path.exists(gsfile):
        x_train_noisy = get_simudata()
        sample_expt, openbeam_expt, transmission_expt = get_data(expt=True)
        sample_phantom, openbeam_phantom, transmission_phantom_gt = get_data(expt=False)
        transmission_phantom = sample_phantom/openbeam_phantom*315715./553690
        transmission_train = x_train_noisy/openbeam_phantom*315715./553690
        g_sample_expt = grad_map_xy(transmission_expt, mask=get_mask(expt=True), agg='p95')
        g_sample_phantom = grad_map_xy(transmission_phantom, mask=get_mask(expt=False), agg='p95')
        g_x_train_noisy = []
        for tidx in np.arange(transmission_train.shape[0]):
            g_x_train_noisy.append(grad_map_xy(transmission_train[tidx], agg='p95'))
        with open(gsfile, 'wb') as f:
            pickle.dump(g_sample_expt, f, 4)
            pickle.dump(g_sample_phantom, f, 4)
            pickle.dump(g_x_train_noisy, f, 4)
    else:
        with open(gsfile, 'rb') as f:
            g_sample_expt = pickle.load(f)
            g_sample_phantom = pickle.load(f)
            g_x_train_noisy = pickle.load(f)

    #g_sample_expt = grad_map_xy(sample_expt)
    #g_sample_phantom = grad_map_xy(sample_phantom)
    #all_vals = np.concatenate([v for v in g_x_train_noisy] + [g_sample_expt, g_sample_phantom])
    all_vals = np.concatenate([g_sample_expt, g_sample_phantom])
    eps = 1e-8
    vmin = np.percentile(all_vals, 0.1)
    vmax = np.percentile(all_vals, 99.9)
    bins = np.exp(np.linspace(np.log(max(vmin, eps)), np.log(vmax), 30))
    c, expt_dens = hist_density(g_sample_expt, bins)
    _, phantom_dens = hist_density(g_sample_phantom,  bins)
    #_, train_dens = hist_density(g_x_train_noisy,  bins)
    train_dens_mat = []
    for v in g_x_train_noisy:
        c, d = hist_density(v, bins)
        train_dens_mat.append(d)
#
    train_dens_mat = np.stack(train_dens_mat, axis=0)
    mean_dens = train_dens_mat.mean(axis=0)
    lo = np.percentile(train_dens_mat, 2.5, axis=0)
    hi = np.percentile(train_dens_mat, 97.5, axis=0)

    plt.figure(figsize=(6,4))
    plt.fill_between(c, lo, hi, color='tab:blue', alpha=0.15, label='Train 95% CI')
    #plt.plot(c, train_dens, color='tab:blue', lw=2, label='Train mean')
    plt.plot(c, mean_dens, color='tab:blue', lw=2, label='Train mean')
    plt.plot(c, phantom_dens, color='tab:orange', lw=2, label='Mock sword')
    plt.plot(c, expt_dens,  color='tab:green',  lw=2, label='Experiment')
    g_train_all = np.concatenate(g_x_train_noisy)
    q95 = np.quantile(g_train_all, 0.95)
    q99 = np.quantile(g_train_all, 0.99)
    plt.axvline(q95, color='r', linestyle='-', lw=1, label='q^train_95', alpha=0.6)
    plt.axvline(q99, color='r', linestyle='--', lw=1, label='q^train_99', alpha=0.6)
    #plt.xscale('log'); 
    #plt.ylim([0, 0.05])
    plt.ylabel('Density'); plt.xlabel('|∇T|')
    plt.legend(); plt.tight_layout(); plt.show()

    grid = np.linspace(q95, np.max([g_sample_phantom.max(), g_sample_expt.max(), g_train_all.max()]), 200)
    x_t, S_t = survival_curve(g_train_all, grid)
    x_m, S_m = survival_curve(g_sample_phantom, grid)
    x_e, S_e = survival_curve(g_sample_expt,  grid)

    plt.figure(figsize=(5.2,3.6))
    plt.semilogy(x_t, S_t, color='tab:blue', label='Train (all)')
    plt.semilogy(x_m, S_m, color='tab:orange', label='Mock sword')
    plt.semilogy(x_e, S_e, color='tab:green',  label='Experiment')
    plt.axvline(q95, color='k', ls='--', alpha=0.5); plt.axvline(q99, color='k', ls=':', alpha=0.5)
    plt.xlabel('|∇T|'); plt.ylabel('Survival  S(g)=P(|∇T|>g)')
    plt.tight_layout(); plt.legend()
    plt.show()


def get_pi_phantom():
    if not os.path.exists('tmp_data_phantom.pkl'):
        tdata, tvar = get_tdata()
        tdata_phantom, tvar_phantom = get_tdata_phantom()
        print(tdata_phantom.shape)
        print(tdata.shape)
        print(tvar_phantom.shape)
        print(tvar.shape)
        M = tdata.shape[0]
        M_phantom = tdata_phantom.shape[0]
        std = (1/M * (tvar.sum(axis=0) + (tdata**2).sum(axis=0)) - np.mean(tdata, axis=0)**2) ** 0.5 
        std_phantom = (1/M_phantom*(tvar_phantom.sum(axis=0) + (tdata_phantom**2).sum(axis=0)) - np.mean(tdata_phantom, axis=0)**2)**0.5
        std_phantom_sample = std_phantom[-2]
        maxstd = np.max(std[:-2], axis=0)
        im = tdata_phantom[0, -2].sum(axis=-1)
        #transmission = tdata[0, -2] / tdata[0, -1]
        #std_sample = std[-2]
        with open('tmp_data_phantom.pkl', 'wb') as f:
            #pickle.dump(transmission, f, 4)
            pickle.dump(im, f, 4)
            pickle.dump(maxstd, f, 4)
            pickle.dump(std_phantom_sample, f, 4)
    else:
        with open('tmp_data_phantom.pkl', 'rb') as f:
            #transmission = pickle.load(f)
            im = pickle.load(f)
            maxstd = pickle.load(f)
            std_phantom_sample = pickle.load(f)
    return std_phantom_sample > maxstd


def get_pi():
    if not os.path.exists('tmp_data.pkl'):
        tdata, tvar = get_tdata()
        M = tdata.shape[0]
        std = (1/M * (tvar.sum(axis=0) + (tdata**2).sum(axis=0)) - np.mean(tdata, axis=0)**2) ** 0.5 
        maxstd = np.max(std[:-2], axis=0)
        im = tdata[0, -2].sum(axis=-1)
        transmission = tdata[0, -2] / tdata[0, -1] 
        std_sample = std[-2]
        with open('tmp_data.pkl', 'wb') as f:
            pickle.dump(transmission, f, 4)
            pickle.dump(im, f, 4)
            pickle.dump(maxstd, f, 4)
            pickle.dump(std_sample, f, 4)
    else:
        with open('tmp_data.pkl', 'rb') as f:
            transmission = pickle.load(f)
            im = pickle.load(f)
            maxstd = pickle.load(f)
            std_sample = pickle.load(f)
    return std_sample>maxstd


def get_gs():
    gsfile = 'gs.pkl'
    if not os.path.exists(gsfile):
        x_train_noisy = get_simudata()
        sample_expt, openbeam_expt, transmission_expt = get_data(expt=True)
        sample_phantom, openbeam_phantom, transmission_phantom_gt = get_data(expt=False)
        transmission_phantom = sample_phantom/openbeam_phantom*315715./553690
        transmission_train = x_train_noisy/openbeam_phantom*315715./553690
        g_sample_expt = grad_map_xy(transmission_expt, mask=get_mask(expt=True), agg='p95')
        g_sample_phantom = grad_map_xy(transmission_phantom, mask=get_mask(expt=False), agg='p95')
        g_x_train_noisy = []
        for tidx in np.arange(transmission_train.shape[0]):
            g_x_train_noisy.append(grad_map_xy(transmission_train[tidx], agg='p95'))
        with open(gsfile, 'wb') as f:
            pickle.dump(g_sample_expt, f, 4)
            pickle.dump(g_sample_phantom, f, 4)
            pickle.dump(g_x_train_noisy, f, 4)
    else:
        with open(gsfile, 'rb') as f:
            g_sample_expt = pickle.load(f)
            g_sample_phantom = pickle.load(f)
            g_x_train_noisy = pickle.load(f)
    return g_sample_expt, g_sample_phantom


def get_g2ds():
    gsfile = 'g2ds.pkl'
    if not os.path.exists(gsfile):
        #x_train_noisy = get_simudata()
        sample_expt, openbeam_expt, transmission_expt = get_data(expt=True)
        sample_phantom, openbeam_phantom, transmission_phantom_gt = get_data(expt=False)
        transmission_phantom = sample_phantom/openbeam_phantom*315715./553690
        #transmission_train = x_train_noisy/openbeam_phantom*315715./553690
        g_sample_expt = grad_map_2d(transmission_expt, mask=get_mask(expt=True), agg='p95')
        g_sample_phantom = grad_map_2d(transmission_phantom, mask=get_mask(expt=False), agg='p95')
        #g_x_train_noisy = []
        #for tidx in np.arange(transmission_train.shape[0]):
        #    g_x_train_noisy.append(grad_map_xy(transmission_train[tidx], agg='p95'))
        with open(gsfile, 'wb') as f:
            pickle.dump(g_sample_expt, f, 4)
            pickle.dump(g_sample_phantom, f, 4)
            #pickle.dump(g_x_train_noisy, f, 4)
    else:
        with open(gsfile, 'rb') as f:
            g_sample_expt = pickle.load(f)
            g_sample_phantom = pickle.load(f)
            #g_x_train_noisy = pickle.load(f)
    return g_sample_expt, g_sample_phantom



#def compare_grads_and_pi():
#    mask = get_mask(expt=True)
#    mask_phantom = get_mask(expt=False)
#    pi_sample = (get_pi().mean(axis=-1))[np.invert(mask)].ravel()
#    pi_sample_phantom = (get_pi_phantom().mean(axis=-1))[np.invert(mask_phantom)].ravel()
#    g_sample, g_sample_phantom = get_gs()




def bin_average(x, y, bins, with_sem=True):
    """
    x: 傾き (flatten)
    y: pi（過小評価発生率, flatten）
    bins: ビン境界（共有）
    with_sem: 標準誤差を返すか
    """
    digit = np.digitize(x, bins)
    xc = 0.5 * (bins[1:] + bins[:-1])  # ビン中心
    y_mean = np.full(len(bins)-1, np.nan)
    y_sem  = np.full(len(bins)-1, np.nan)

    for i in range(1, len(bins)):
        sel = (digit == i)
        if np.any(sel):
            vals = y[sel]
            y_mean[i-1] = np.nanmean(vals)
            if with_sem:
                # 標準誤差 (n有効のとき)
                n = np.sum(~np.isnan(vals))
                if n > 1:
                    y_sem[i-1] = np.nanstd(vals, ddof=1) / np.sqrt(n)
                elif n == 1:
                    y_sem[i-1] = 0.0
    return xc, y_mean, y_sem




from scipy.ndimage import gaussian_filter, map_coordinates, distance_transform_edt
from skimage.morphology import medial_axis
from skimage.measure import label, regionprops

# ------------------------------------------------------------
# 1) 棟近傍の「中心線」（おおよその ridge）を抽出
#    - 入力 mask(True=背景) から有効領域(~mask)の medial axis をとる
#    - g2d の強い場所（上位 q%）の近傍にある軸点だけを残して「棟の中心線」を近似
# ------------------------------------------------------------
def extract_centerline(mask, g2d, ridge_q=0.80, smooth_sigma=1.0, min_points=100):
    """
    mask: True=background
    g2d : |∇T| 2D
    ridge_q: g2d の上位分位（例: 0.8）を「崖領域」とみなし、その近傍にある軸点を優先
    smooth_sigma: g2d を少し平滑化して安定化
    return:
        pts (N,2) ndarray in (y,x) sorted along principal axis
    """
    H, W = g2d.shape
    g2 = gaussian_filter(g2d, smooth_sigma) if (smooth_sigma and smooth_sigma>0) else g2d.copy()
    valid = ~mask

    # medial axis（中心線）と距離画像
    skel, dist = medial_axis(valid, return_distance=True)

    # 崖（高 g）領域（上位分位）
    th = np.quantile(g2[valid], ridge_q)
    ridge_zone = np.zeros_like(g2, dtype=bool)
    ridge_zone[valid] = g2[valid] >= th

    # 中心線の中でも「崖領域に近い」点のみ抽出（膨張的に許容）
    # ここは EDT を使って ridge_zone からの距離が小さい軸点を選ぶ
    d2ridge = distance_transform_edt(~ridge_zone)
    sel = skel & (d2ridge <= 2.0)  # 2px 以内（適宜調整）

    # 保険：もし選ばれた点が極端に少なければ skel そのものを採用
    if sel.sum() < min_points:
        sel = skel.copy()

    ys, xs = np.nonzero(sel)
    pts = np.stack([ys, xs], axis=1)
    if len(pts) == 0:
        # 中心線が取れなかった場合：有効領域の最長主軸をざっくり返す
        ys, xs = np.nonzero(valid)
        pts = np.stack([ys, xs], axis=1)

    # 主成分方向で並べ替えて「棟の長手方向に沿った順序」を与える
    # （2D PCA：符号任意性はあるが、ソートには影響しない）
    P = pts.astype(float)
    c = P.mean(axis=0)
    X = P - c
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    v = Vt[0]  # 第一主成分方向（2次元）
    # 射影値でソート
    t = X @ v
    order = np.argsort(t)
    pts_sorted = pts[order]
    return pts_sorted  # (N,2) in (y,x)


def extract_ridge_from_mask_simple(mask, take='top', min_step=2):
    """
    mask: True=背景 / False=刀本体
    take: 'top'（上側）/ 'bottom'（下側）から選ぶだけの簡易版
    min_step: 近接点の間引き（px）
    return: pts (N,2) in (y,x) sorted along x 昇順
    """
    valid = ~mask  # 刀本体
    # 1px 幅の外形（XOR で境界線）
    boundary = binary_dilation(valid) ^ valid

    H, W = boundary.shape
    ys, xs = [], []

    # 列ごとに“最上位”または“最下位”の境界画素を拾う
    for x in range(W):
        col = boundary[:, x].nonzero()[0]
        if col.size == 0:
            continue
        y = col.min() if take == 'top' else col.max()
        ys.append(y); xs.append(x)

    pts = np.stack([np.array(ys), np.array(xs)], axis=1) if len(ys) > 0 else np.zeros((0,2), int)

    # 近接点を軽く間引き（ギザギザ抑制）
    keep = [True]
    for i in range(1, len(pts)):
        if np.hypot(*(pts[i] - pts[i-1])) < min_step:
            keep.append(False)
        else:
            keep.append(True)
    if len(pts) > 0:
        pts = pts[np.array(keep)]

    return pts


def extract_spine_from_mask_by_edt(mask, min_step=2, smooth_sigma=1.0):
    """
    mask: True=背景 / False=刀本体
    return: ctr (N,2) in (y,x) : 棟（spine）に相当する境界の点列（xで単調）
    ロジック:
      1) valid = ~mask（刀領域）から 1px 幅の外形 boundary を抽出
      2) 列(x)ごとに「最上位(top)」と「最下位(bottom)」の境界画素を拾う
      3) valid の距離変換 EDT を計算し、各候補境界の EDT を比較
         → 平均/中央値が大きい側を「棟」とみなす
    """
    valid = ~mask
    boundary = binary_dilation(valid) ^ valid
    H, W = boundary.shape

    # 列ごとに top/bottom の境界点を拾う
    top_pts, bot_pts = [], []
    for x in range(W):
        ys = np.nonzero(boundary[:, x])[0]
        if ys.size == 0:
            continue
        top_pts.append((ys.min(), x))
        bot_pts.append((ys.max(), x))
    if len(top_pts) == 0 or len(bot_pts) == 0:
        return np.zeros((0,2), dtype=int)

    top_pts = np.array(top_pts, dtype=int)
    bot_pts = np.array(bot_pts, dtype=int)

    # EDT（内部からの距離）。内部(valid)のみ >0
    edt = distance_transform_edt(valid)

    # それぞれの境界点での EDT 値を取得（列ごと対応）
    top_edt = edt[top_pts[:,0], top_pts[:,1]]
    bot_edt = edt[bot_pts[:,0], bot_pts[:,1]]

    # 「棟＝厚い側」→ EDT が大きい方を採用（中央値で判定が安定）
    if np.nanmedian(top_edt) >= np.nanmedian(bot_edt):
        ctr = top_pts.copy()
    else:
        ctr = bot_pts.copy()

    # 近接間引き（ギザ抑制）
    keep = [True]
    for i in range(1, len(ctr)):
        if np.hypot(*(ctr[i] - ctr[i-1])) < min_step:
            keep.append(False)
        else:
            keep.append(True)
    ctr = ctr[np.array(keep)]

    # y を軽く平滑（オプション）
    if smooth_sigma and len(ctr) > 3:
        ys = ctr[:,0].astype(float)
        xs = ctr[:,1].astype(float)
        ys_sm = gaussian_filter1d(ys, sigma=smooth_sigma)
        ctr = np.stack([ys_sm, xs], axis=1)

    return ctr.astype(int)



# ------------------------------------------------------------
# 2) 中心線に沿った接線→法線を計算し、法線方向に g2d をサンプリング
# ------------------------------------------------------------
def sample_profiles_along_normals(g2d, centerline_pts, t_half_width=20.0, dt=0.5, 
                                  keep_fraction=1.0, min_sep_px=3):
    """
    g2d: 2D
    centerline_pts: (N,2) (y,x) sorted
    t_half_width: 片側サンプリング距離（px）
    dt: サンプリング間隔（px）
    keep_fraction: 中心線点の間引き率（0< f <=1, 例:0.5なら半分だけ使う）
    min_sep_px: 近接しすぎる点を間引く最小距離
    return:
        profiles: (M, T)  1点あたりの法線方向サンプル
        t_axis  : (T,)    法線距離（負/正を含む）
    """
    H, W = g2d.shape
    pts = centerline_pts.copy()

    # 近接点の間引き（粗く）
    keep = [True]
    for i in range(1, len(pts)):
        if np.hypot(*(pts[i] - pts[i-1])) < min_sep_px:
            keep.append(False)
        else:
            keep.append(True)
    pts = pts[np.array(keep)]

    ## ランダムにさらに間引き
    #if keep_fraction < 1.0:
    #    rng = np.random.default_rng(42)
    #    idx = rng.choice(len(pts), size=max(5, int(len(pts)*keep_fraction)), replace=False)
    #    idx.sort()
    #    pts = pts[idx]

    # --- ランダム間引き（安全版） ---
    n = len(pts)
    if keep_fraction < 1.0 and n > 0:
        target = int(np.ceil(n * keep_fraction))
        # 常に 1..n にクリップ
        target = max(1, min(n, target))
        if target < n:
            rng = np.random.default_rng(42)
            idx = rng.choice(n, size=target, replace=False)
            idx.sort()
            pts = pts[idx]
# ----------------------------------


    # 接線ベクトル（前後差分）→ 法線（90度回転）
#    def tangent(i):
#        if i == 0:
#            v = pts[i+1] - pts[i]
#        elif i == len(pts)-1:
#            v = pts[i] - pts[i-1]
#        else:
#            v = pts[i+1] - pts[i-1]
#        v = v.astype(float)
#        nrm = np.linalg.norm(v)
#        if nrm == 0: 
#            return np.array([0.0, 1.0])
#        return v / nrm


    # --- inside sample_profiles_along_normals() ---
    def tangent(i):
        n = len(pts)
        if n == 1:
            # 1点しかない場合のフォールバック（単位ベクトルを返す）
            return np.array([0.0, 1.0], dtype=float)

        if i == 0:
            v = pts[i+1] - pts[i]
        elif i == n - 1:
            v = pts[i] - pts[i-1]
        else:
            v = pts[i+1] - pts[i-1]
        v = v.astype(float)
        nrm = np.linalg.norm(v)
        if nrm == 0:
            return np.array([0.0, 1.0], dtype=float)
        return v / nrm


    # 法線
    def normal_from_tan(tvec):
        # 2D で (ty, tx) に対し法線の一つは (-tx, ty)
        return np.array([-tvec[1], tvec[0]])

    t_axis = np.arange(-t_half_width, t_half_width + 1e-6, dt)  # 含む
    profiles = []
    for i, (yy, xx) in enumerate(pts):
        tvec = tangent(i)
        nvec = normal_from_tan(tvec)
        # 法線方向座標列
        ys = yy + t_axis * nvec[0]
        xs = xx + t_axis * nvec[1]
        # 画像外を避けて補間
        coords = np.vstack([ys, xs])
        prof = map_coordinates(g2d, coords, order=1, mode='nearest')
        profiles.append(prof)
    profiles = np.asarray(profiles)  # (M, T)
    return profiles, t_axis

# ------------------------------------------------------------
# 3) 平均プロファイルと SEM の計算
# ------------------------------------------------------------
def mean_and_sem(profiles, axis=0):
    """
    profiles: (M, T) など
    axis: 平均化する軸（通常0 = 本数方向）
    return: mean, sem
    """
    mean = np.nanmean(profiles, axis=axis)
    sd   = np.nanstd(profiles,  axis=axis, ddof=1)
    n    = np.sum(~np.isnan(profiles), axis=axis)
    sem  = sd / np.sqrt(np.maximum(n, 1))
    return mean, sem

# ------------------------------------------------------------
# 4) 上記をまとめて 実験/模擬刀 の平均±SEM プロファイルを描画
# ------------------------------------------------------------
def plot_ridge_normal_profiles_sem(g2d_exp, g2d_ph, mask_exp=None, mask_ph=None,
                                   ridge_q=0.80, t_half_width=20.0, dt=0.5,
                                   keep_fraction=0.6, min_sep_px=3,
                                   smooth_sigma_centerline=1.0,
                                   smooth_sigma_signal=0.5,
                                   use_ci95=False, title='Ridge-normal profile (mean ± SEM)'):
    """
    use_ci95: True にすると ±1.96*SEM（95%CI相当）をバンドにする。False なら ±SEM。
    """
    if mask_exp is None:
        mask_exp = np.zeros_like(g2d_exp, dtype=bool)
    if mask_ph is None:
        mask_ph = np.zeros_like(g2d_ph, dtype=bool)

    # 1) centerline
    #ctr_exp = extract_centerline(mask_exp, g2d_exp, ridge_q=ridge_q, smooth_sigma=smooth_sigma_centerline)
    #ctr_ph  = extract_centerline(mask_ph,  g2d_ph,  ridge_q=ridge_q, smooth_sigma=smooth_sigma_centerline)

    # 1) centerline を mask の外形境界（棟側）から単純に取得
    #ctr_exp = extract_ridge_from_mask_simple(mask_exp, take='top', min_step=min_sep_px)
    #ctr_ph  = extract_ridge_from_mask_simple(mask_ph,  take='top', min_step=min_sep_px)


    # ↓ これに置き換え
    ctr_exp = extract_spine_from_mask_by_edt(mask_exp, min_step=min_sep_px, smooth_sigma=1.0)
    ctr_ph  = extract_spine_from_mask_by_edt(mask_ph,  min_step=min_sep_px, smooth_sigma=1.0)


    # 2) sampling
    profE, t_axis = sample_profiles_along_normals(
        gaussian_filter(g2d_exp, smooth_sigma_signal) if smooth_sigma_signal>0 else g2d_exp,
        ctr_exp, t_half_width=t_half_width, dt=dt, keep_fraction=keep_fraction, min_sep_px=min_sep_px
    )
    profP, _ = sample_profiles_along_normals(
        gaussian_filter(g2d_ph, smooth_sigma_signal) if smooth_sigma_signal>0 else g2d_ph,
        ctr_ph, t_half_width=t_half_width, dt=dt, keep_fraction=keep_fraction, min_sep_px=min_sep_px
    )

    # 3) mean & SEM
    meanE, semE = mean_and_sem(profiles=profE, axis=0)
    meanP, semP = mean_and_sem(profiles=profP, axis=0)

    k = 1.96 if use_ci95 else 1.0
    loE, hiE = meanE - k*semE, meanE + k*semE
    loP, hiP = meanP - k*semP, meanP + k*semP

    # 4) plot
    plt.figure(figsize=(7,5))
    plt.plot(t_axis, meanE, color='tab:blue',  label='Experiment')
    plt.fill_between(t_axis, loE, hiE, color='tab:blue', alpha=0.25, linewidth=0)

    plt.plot(t_axis, meanP, color='tab:orange', label='Phantom')
    plt.fill_between(t_axis, loP, hiP, color='tab:orange', alpha=0.25, linewidth=0)

    plt.axvline(0, color='k', lw=0.8, alpha=0.3)
    plt.xlabel('Normal distance t [px]')
    ylabel = r'$|\nabla T(x,y)|$ (g2d)'
    if use_ci95:
        ylabel += '  (mean ± 95%CI)'
    else:
        ylabel += '  (mean ± SEM)'
    plt.ylabel(ylabel)
    plt.title(title + ('' if not use_ci95 else ' — 95%CI'))
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()


    fig, ax = plt.subplots(2)
    ax[0].imshow(g2d_ph)
    ax[0].plot(ctr_ph[:,1], ctr_ph[:,0], c='r')
    ax[1].imshow(g2d_exp)
    ax[1].plot(ctr_exp[:,1], ctr_exp[:,0], c='r')
    plt.show()


    # 返り値（必要なら解析に）
    return {
        't_axis': t_axis,
        'exp': {'profiles': profE, 'mean': meanE, 'sem': semE},
        'ph' : {'profiles': profP, 'mean': meanP, 'sem': semP},
        'centerline': {'exp': ctr_exp, 'ph': ctr_ph}
    }




import numpy as np
from scipy.ndimage import map_coordinates

def _first_exit_along_ray(valid, y0, x0, vy, vx, step=0.5, max_steps=2000):
    """
    valid: True=刀内部（~mask）
    (y0,x0): 出発点（中心線上）
    (vy,vx): 単位ベクトル（進行方向）
    step: ステップ長 [px]
    return: (y_end, x_end, traveled_length)
      - 出発点を含む内部から出て最初に「外部 or 画外」になる直前の点を返す
      - 内部がすぐ途切れる場合は出発点を返す
    """
    H, W = valid.shape
    y, x = float(y0), float(x0)
    L = 0.0
    for k in range(int(max_steps)):
        yn = y + step * vy
        xn = x + step * vx
        # 画外に出たら終了（直前位置を端点とみなす）
        if not (0 <= yn < H and 0 <= xn < W):
            return y, x, L
        # 次の位置が内部かチェック
        if not valid[int(round(yn)), int(round(xn))]:
            # 外に出る直前の位置を端点とみなす
            return y, x, L
        # 進行
        y, x = yn, xn
        L += step
    return y, x, L


def sample_profiles_full_crossing(g2d, mask, centerline_pts, dt=0.5, step_edge=0.5,
                                  normalize_to_unit=True, min_span_px=5.0):
    """
    g2d: 2D ndarray
    mask: True=背景 / False=刀本体
    centerline_pts: (N,2) in (y,x) - 既に得られた中心線（棟）上の点列（順序は任意）
    dt: 線分上サンプル間隔 [px]
    step_edge: 端点探索のステップ長 [px]
    normalize_to_unit: Trueなら各プロファイルを [-1, +1] の共通座標に再サンプリング
    min_span_px: 刃～棟の物理スパンがこの値未満なら破棄（ノイズ除去）

    return:
      profiles_px  : list of 1D arrays（各プロファイルの物理ピクセル座標での g2d 値）
      t_axes_px    : list of 1D arrays（各プロファイル固有の t 軸 [px]、t=0 は中心線位置）
      profiles_unit: 2D ndarray (M, T)  もし normalize_to_unit=True の場合、
                     すべてを [-1, +1] 共通座標に再サンプル済み（M=本数、T=共通点数）
      t_axis_unit  : 1D ndarray 共通の単位軸 [-1..+1]
      meta         : dict（各プロファイルの端点座標など）
    """
    H, W = g2d.shape
    valid = ~mask  # True=内部

    # 中心線点を x（列）でソートしておくと描画が綺麗
    # （順序に意味がない場合は無視してよい）
    pts = np.array(centerline_pts, dtype=float)
    if len(pts) == 0:
        return [], [], None, None, {}

    # 進行ごとに結果を格納
    profiles_px = []
    t_axes_px = []
    endpoints_a = []  # 負方向端点
    endpoints_b = []  # 正方向端点
    spans = []        # 端点間の全長 [px]

    for i, (yy, xx) in enumerate(pts):
        # 局所接線（前後差分）
        if i == 0:
            v = pts[i+1] - pts[i]
        elif i == len(pts)-1:
            v = pts[i] - pts[i-1]
        else:
            v = pts[i+1] - pts[i-1]
        v = v.astype(float)
        nrm = np.linalg.norm(v)
        if nrm == 0:
            # デフォルト向き
            tvec = np.array([0.0, 1.0])
        else:
            tvec = v / nrm

        # 法線（90°回転）
        nvec = np.array([-tvec[1], tvec[0]], dtype=float)

        # 法線の両方向に「内部→外部」の境界まで伸ばす
        ya, xa, La = _first_exit_along_ray(valid, yy, xx, -nvec[0], -nvec[1], step=step_edge)
        yb, xb, Lb = _first_exit_along_ray(valid, yy, xx,  nvec[0],  nvec[1], step=step_edge)

        span = La + Lb
        if span < min_span_px:
            continue  # 刃～棟距離が短すぎる点はスキップ（先端など）

        # 端点から端点までサンプル（t=0 を中心線位置にし、[-La, +Lb] を等間隔）
        t_axis = np.arange(-La, Lb + 1e-6, dt)
        ys = yy + t_axis * nvec[0]
        xs = xx + t_axis * nvec[1]
        coords = np.vstack([ys, xs])
        prof = map_coordinates(g2d, coords, order=1, mode='nearest')

        profiles_px.append(prof)
        t_axes_px.append(t_axis)
        endpoints_a.append((ya, xa))
        endpoints_b.append((yb, xb))
        spans.append(span)

    meta = {
        "endpoints_neg": np.array(endpoints_a, dtype=float),
        "endpoints_pos": np.array(endpoints_b, dtype=float),
        "spans_px": np.array(spans, dtype=float),
        "centerline": pts
    }

    # ここから共通座標（-1..+1）への正規化（任意）
    profiles_unit = None
    t_axis_unit = None
    if len(profiles_px) > 0 and normalize_to_unit:
        # 各プロファイルごとに独自の [-La, +Lb] を [-1, +1] へ線形写像して再サンプリング
        T = 201  # 共通グリッドの点数（適宜）
        t_axis_unit = np.linspace(-1.0, 1.0, T)
        profiles_u = []
        for prof, tpx in zip(profiles_px, t_axes_px):
            # 物理→単位：u = (t - tmin)/(tmax - tmin) * 2 - 1
            tmin, tmax = tpx[0], tpx[-1]
            u_of_t = (tpx - tmin) / (tmax - tmin) * 2.0 - 1.0
            # u の単調性は保証されるので 1D 線形補間
            prof_u = np.interp(t_axis_unit, u_of_t, prof, left=np.nan, right=np.nan)
            profiles_u.append(prof_u)
        profiles_unit = np.array(profiles_u)

    return profiles_px, t_axes_px, profiles_unit, t_axis_unit, meta





def densify_centerline_by_columns(ctr, width, x_margin=0):
    """
    ctr: (N,2) in (y,x) の中心線点列
    width: 画像の幅 W
    x_margin: 左右マージン（外形の端の不安定列を除外したいときに>0）
    return: ctr_dense (M,2) in (y,x) - 列ごとに 1 点（x は整数）
    """
    ctr = np.asarray(ctr, dtype=float)
    if ctr.size == 0:
        return np.zeros((0,2), float)

    # xでソート
    order = np.argsort(ctr[:,1])
    xs = ctr[order, 1]
    ys = ctr[order, 0]

    # x の範囲
    x_min = max(int(np.ceil(xs.min())), x_margin)
    x_max = min(int(np.floor(xs.max())), width - 1 - x_margin)
    if x_max <= x_min:
        return np.zeros((0,2), float)

    # 欠損列も含めて y(x) を線形補間（単調でなくてOK）
    x_grid = np.arange(x_min, x_max + 1)  # すべての整数列
    y_grid = np.interp(x_grid, xs, ys)    # 1D 線形補間

    ctr_dense = np.stack([y_grid, x_grid], axis=1)
    return ctr_dense


def compare_grads_and_pi():

    mask = get_mask(expt=True)
    mask_phantom = get_mask(expt=False)

    # 1) 1D 化された pi
    pi_sample = (get_pi().mean(axis=-1))[~mask].ravel()
    pi_sample_phantom = (get_pi_phantom().mean(axis=-1))[~mask_phantom].ravel()

    # 2) 1D 化された g（前処理側で mask 済みだと仮定。未適用なら同様に ~mask で抽出）
    g_sample, g_sample_phantom = get_gs()  # flatten 済みが返る想定
    g2d_sample, g2d_sample_phantom = get_g2ds()

    fig, ax = plt.subplots(3)
    import cv2
    width = g2d_sample.shape[0]
    height = g2d_sample.shape[1]
    center = (int(width/2), int(height/2))
    angle = 5.5
    angle_ph = -7.
    scale = 1.
    trans = cv2.getRotationMatrix2D(center, angle, scale)
    trans_ph =  cv2.getRotationMatrix2D(center, angle_ph, scale)
    transp1 = cv2.getRotationMatrix2D(center, angle+0.5, scale)
    transm1 = cv2.getRotationMatrix2D(center, angle-0.5, scale)
    trans_php1 = cv2.getRotationMatrix2D(center, angle_ph+0.25, scale)
    trans_phm1 = cv2.getRotationMatrix2D(center, angle_ph-0.25, scale)
    g2d_sample_rotated = cv2.warpAffine(g2d_sample, trans, (height, width))
    g2d_sample_p1 = cv2.warpAffine(g2d_sample, transp1, (height, width))
    g2d_sample_m1 = cv2.warpAffine(g2d_sample, transm1, (height, width))
    g2d_sample_php1 = cv2.warpAffine(g2d_sample_phantom, trans_php1, (height, width))[::-1, :]
    g2d_sample_phm1 = cv2.warpAffine(g2d_sample_phantom, trans_phm1, (height, width))[::-1, :]
    g2d_sample_phantom_rotated = cv2.warpAffine(g2d_sample_phantom, trans_ph, (height, width))[::-1, :]
    ax[2].plot(g2d_sample_rotated.mean(axis=-1), label='expt')
    #plt.plot(g2d_sample_php1.mean(axis=-1), label='phantomp1', lw=0.5)
    #plt.plot(g2d_sample_phm1.mean(axis=-1), label='phantomm1', lw=0.5)
    ax[2].plot(g2d_sample_phantom_rotated.mean(axis=-1), label='phantom')
    ax[2].set_xlim([0, 50])

    vmin = min(np.mean(g2d_sample), np.mean(g2d_sample_phantom))
    vmax = max(np.max(g2d_sample), np.max(g2d_sample_phantom))
    
    im1=ax[0].imshow(g2d_sample_rotated, vmin=vmin, vmax=vmax)
    im2=ax[1].imshow(g2d_sample_phantom_rotated, vmin=vmin, vmax=vmax)
    plt.legend()
    plt.show()

    cbar = fig.colorbar(im1, ax=ax.ravel().tolist(), fraction=0.025, pad=0.02)
    cbar.set_label(r'$|\nabla T(x,y)|$ (Sobel)')

    #plt.colorbar(im)
    #plt.colorbar(im2)
    plt.show()


    # 3) 共有ビン（両方の範囲をカバー）
    g_min = np.nanmin([np.nanmin(g_sample), np.nanmin(g_sample_phantom)])
    g_max = np.nanmax([np.nanmax(g_sample), np.nanmax(g_sample_phantom)])
    bins = np.linspace(g_min, g_max, 41)  # 40ビン（必要に応じて変更）

    # 4) ビン平均（平均と標準誤差）
    xc_s, mean_s, sem_s = bin_average(g_sample, pi_sample, bins)
    xc_p, mean_p, sem_p = bin_average(g_sample_phantom, pi_sample_phantom, bins)

    # 5) プロット（論文用の体裁）
    plt.figure(figsize=(7,5))
    plt.errorbar(xc_s, mean_s, yerr=sem_s, fmt='-o', label='Experiment', capsize=3)
    plt.errorbar(xc_p, mean_p, yerr=sem_p, fmt='-s', label='Phantom', capsize=3)
    plt.xlabel('Transmission gradient |∇T| (Sobel, smoothed)')
    plt.ylabel('Underestimation rate π')
    plt.title('Correlation between gradient and underestimation (π)')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig('fig_bi_g_vs_pi.png')
    plt.show()
    return bins, (xc_s, mean_s, sem_s), (xc_p, mean_p, sem_p)






if __name__ == '__main__':
    compare_grads()
    compare_grads_and_pi()

