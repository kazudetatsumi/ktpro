#!?usr/bin/env python
import pickle
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.gridspec import GridSpec
import numpy as np
import os
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
    import copy
    _g2d = copy.deepcopy(g2d)
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
    plt.ylim([0, 0.05])
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



if __name__ == '__main__':
    compare_grads()


