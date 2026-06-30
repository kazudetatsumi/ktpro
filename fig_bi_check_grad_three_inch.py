#!/usr/bin/env python
# integrate_grad_onefig_inches.py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import socket
import os

# use exsiting preprocess and data
import fig_bi_check_grad as mod

# setting on fonts and eps file
mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 9/2*1.3, "axes.labelsize": 11/2*1.3, "xtick.labelsize": 9/2*1.3, "ytick.labelsize": 9/2*1.3,
    "legend.fontsize": 9/2*1.3,
    'xtick.major.pad': 2.2,
    'ytick.major.pad': 2.2,
})
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'


def draw_density(ax, gsfile):
    # histogram of transmission gradient
    # shade and mean are drawn if gsfile contains whole train data points
    g_expt, g_ph = mod.get_gs(gsfile=gsfile)
    all_vals = np.concatenate([g_expt, g_ph])
    eps = 1e-8
    vmin = np.percentile(all_vals, 0.1)
    vmax = np.percentile(all_vals, 99.9)
    bins = np.exp(np.linspace(np.log(max(vmin, eps)), np.log(vmax), 100))

    c, expt_dens = mod.hist_density(g_expt, bins)
    _, ph_dens = mod.hist_density(g_ph, bins)

    # 95%CI and mean are used if whole train data points available.
    mean_dens = lo = hi = None
    try:
        import pickle, os
        if os.path.exists(gsfile):
            with open(gsfile, 'rb') as f:
                _ = pickle.load(f); _ = pickle.load(f); g_x_train_noisy = pickle.load(f)
            mat = []
            for v in g_x_train_noisy:
                _, d = mod.hist_density(v, bins)
                mat.append(d)
            mat = np.stack(mat, 0)
            mean_dens = mat.mean(0)
            lo = np.percentile(mat, 2.5,  axis=0)
            hi = np.percentile(mat, 97.5, axis=0)
    except Exception:
        pass

    if lo is not None and hi is not None:
        ax.fill_between(c, lo, hi, color='lightgray', alpha=0.45,
                        label='Train CI')
    if mean_dens is not None:
        ax.plot(c, mean_dens, color='gray', lw=1.5, label='Train')

    ax.plot(c, ph_dens,   lw=1.5, label='Moc', ls='--', c='k')
    ax.plot(c, expt_dens, lw=1.5, label='Expt', ls='-', c='k')
    # ax.set_ylim(0, 0.05)
    ax.set_yscale('log')
    ax.set_xlabel(r'$|\nabla T|_{qtof95}$')
    ax.set_ylabel('Probability density')
    # ax.legend(loc='best', frameon=False)
    ax.legend(bbox_to_anchor=(1.01, 1.03), loc='upper right',
              ncol=2, frameon=False, borderpad=0.1, labelspacing=0.1,
              handlelength=2., handletextpad=0.5, columnspacing=0.5,)
    # ax.tick_params(axis='both', which='both', pad=0.4)
    #ax.set_title(r'Probability Density of $|\nabla T|_{qtof95}$')


def draw_corr(
        ax,
        gsfile,
        tmp_datafile,
        tmp_dataphantomfile,
        tdatafile,
        tdataphantomfile,
        ):
    # gradient vs violation prob correlation
    mask = mod.get_mask(expt=True)
    mask_ph = mod.get_mask(expt=False)
    pi_s = (mod.get_pi(
               tmp_datafile=tmp_datafile,
               tdatafile=tdatafile,
                       ).mean(axis=-1)
            )[~mask].ravel()
    pi_ph = (mod.get_pi_phantom(
                tmp_dataphantomfile=tmp_dataphantomfile,
                tdatafile=tdatafile,
                tdataphantomfile=tdataphantomfile,
                                ).mean(axis=-1)
             )[~mask_ph].ravel()
    g_s, g_ph = mod.get_gs(gsfile=gsfile)

    g_min = np.nanmin([np.nanmin(g_s),  np.nanmin(g_ph)])
    g_max = np.nanmax([np.nanmax(g_s),  np.nanmax(g_ph)])
    bins = np.linspace(g_min, g_max, 41)

    xc_s, mean_s, sem_s = mod.bin_average(g_s,  pi_s,  bins)
    xc_p, mean_p, sem_p = mod.bin_average(g_ph, pi_ph, bins)

    ax.errorbar(xc_s, mean_s, yerr=sem_s, color='black', marker='o', lw=0.7,
                ls='-', label='Experiment', capsize=2, ms=2.5)
    ax.errorbar(xc_p, mean_p, yerr=sem_p, color='black', marker='s', lw=0.7,
                ls='-', mfc='white', mew=1, label='Phantom', capsize=2, ms=2.5)
    ax.set_xlabel(r'Transmission gradient $|\nabla T|_{qtof95}$', labelpad=-0.5)
    ax.set_ylabel('Violation rate')
    #ax.set_title(r'Correlation between $|\nabla T|_{qtof95}$ and π')
    ax.legend(loc='best', frameon=False)


def draw_line(ax, gdsfile):
    # Draw row mean line plots on spatially rotated map
    import cv2
    g2d_s, g2d_ph = mod.get_g2ds(gdsfile=gdsfile)
    H, W = g2d_s.shape
    center = (int(H/2), int(W/2))
    angle, angle_ph = 5.5, -7.0
    scale = 1.0
    trans = cv2.getRotationMatrix2D(center, angle,    scale)
    trans_ph = cv2.getRotationMatrix2D(center, angle_ph, scale)
    g2d_s_rot = cv2.warpAffine(g2d_s,  trans,    (W, H))
    g2d_ph_rot = cv2.warpAffine(g2d_ph, trans_ph, (W, H))[::-1, :]

    ax.plot(g2d_s_rot.mean(axis=-1), ls='-', label='Experiment', c='k')
    ax.plot(g2d_ph_rot.mean(axis=-1), ls='--', label='Phantom', c='k')
    ax.set_xlim(0, 50)
    ax.set_xlabel('Transverse position / ch', labelpad=-0.5)
    ax.set_ylabel(r'$<|\nabla T|_{qtof95}>_{longitude}$')
    ax.legend(loc='best', frameon=False)
    #ax.set_title(r'Spatial Distribution of $<|\nabla T|_{qtof95}>_{longitude}$')


def get_axes_and_fig():
    fig = plt.figure()
    # strict inch scale
    left = 0.65/2
    right = 0.45/2
    top = 0.35/2
    bottom = 0.55/2
    panel_w = 4.5/2.5
    h1 = 1.55/1.9
    h2 = 1.55/1.9
    h3 = 1.55/1.9
    # h4 = 1.25
    s12 = 0.28*2.5/2.
    s23 = 0.25*2.5/2.
    # s34 = 0.28*2
    """
    1列×4段の Axes をインチ指定で厳密配置し、[ax1, ax2, ax3, ax4] を返す。
    """
    # 図全体サイズ（inch）
    w = left + panel_w + right
    #h = top + (h1 + s12 + h2 + s23 + h3 + s34 + h4) + bottom
    h = top + (h1 + s12 + h2 + s23 + h3) + bottom
    # 図の物理サイズを決め打ち
    print(w, h)
    fig.set_size_inches(w, h)

    # 下から段を積み上げる
    _y = bottom
    # パネル4（最下段）
    x0 = left/w
    w = panel_w/w
    ax3 = fig.add_axes([x0, _y/h, w, h3/h])
    _y += h3 + s23
    # パネル3
    ax2 = fig.add_axes([x0, _y/h, w, h2/h])
    _y += h2 + s12
    # パネル2
    ax1 = fig.add_axes([x0, _y/h, w, h1/h])
    #_y += h2 + s12
    # パネル1（最上段）
    #ax1 = fig.add_axes([x0, _y/h, w, h1/h])
    # ax1, ax2, ax3, ax4 = axs  # ax1=最上段 … ax4=最下段
    return fig, ax1, ax2, ax3


def main():
    hostname = socket.gethostname()
    if hostname == "mlfdev61":
        fdir = "/home/kazu/restormer_rev2_lim/bi3d/restormer_conv3d/" +\
             "for_single/train/full/211/true_edge/nll/gau2ch/ktrand/rev4/"
    elif hostname == "mlfdev51":
        fdir = "/data1/kazu/restormer_rev2_lim/bi3d/restormer_conv3d/" +\
             "for_single/train/full/211/true_edge/nll/gau2ch/ktrand/"
    else:
        raise RuntimeError("You should be in mlfdev51 or mlfdev61")
    gsfile = fdir + 'gs.pkl'
    if not os.path.exists(gsfile):
        raise RuntimeError("gsfile doest not exist")
    gdsfile = fdir + 'g2ds.pkl'
    if not os.path.exists(gdsfile):
        raise RuntimeError("gdsfile doest not exist")
    tmp_datafile = fdir + 'tmp_data.pkl'
    if not os.path.exists(tmp_datafile):
        raise RuntimeError("tmp_datafile doest not exist")
    tmp_dataphantomfile = fdir + 'tmp_data_phantom.pkl'
    if not os.path.exists(tmp_dataphantomfile):
        raise RuntimeError("tmp_dataphantomfile doest not exist")
    tdatafile = fdir + 'tdata.pkl'
    if (not os.path.exists(tdatafile)) and (not os.path.exists(tmp_datafile)):
        raise RuntimeError("tdatafile and tmp_datafile do not exist")
    tdataphantomfile = fdir + 'tdata_phantom.pkl'
    if not (os.path.exists(tdataphantomfile)) and\
           (not os.path.exists(tmp_dataphantomfile)):
        raise RuntimeError(tdataphantomfile + " and " + tmp_dataphantomfile +
                           " do not exist")
    fig, ax1, ax2, ax3 = get_axes_and_fig()
    # 並べ順は (1) → (2) → (4) → (3)
    draw_density(ax2, gsfile)
    #draw_survival(ax3)
    draw_corr(
            ax1,
            gsfile,
            tmp_datafile,
            tmp_dataphantomfile,
            tdatafile,
            tdataphantomfile,
            )   # ← 元「4番目」
    draw_line(ax3, gdsfile)   # ← 元「3番目（マップなし）」

    # 体裁：上3段の x ラベルは隠す（横位置は inch 固定なのでズレません）
    #for a in [ax1, ax2, ax3]:
    #    a.tick_params(labelbottom=False)

    fig.savefig("fig_bi_grad.eps")
    plt.show()
    #print(f"Saved -> {output}")


if __name__ == "__main__":
    main()
