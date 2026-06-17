#!/usr/bin/env python
# integrate_grad_onefig_inches.py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# 既存の計算・前処理を流用
import fig_bi_check_grad as mod   # ← あなたのファイル名

# フォント・EPS設定
mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 9, "axes.labelsize": 11, "xtick.labelsize": 9, "ytick.labelsize": 9,
    "legend.fontsize": 9,
})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'


def draw_density(ax):
    # 勾配ヒスト（Train帯が gs.pkl にあればシェードと平均も描く）
    g_expt, g_ph = mod.get_gs()
    all_vals = np.concatenate([g_expt, g_ph])
    eps = 1e-8
    vmin = np.percentile(all_vals, 0.1)
    vmax = np.percentile(all_vals, 99.9)
    bins = np.exp(np.linspace(np.log(max(vmin, eps)), np.log(vmax), 100))

    c, expt_dens = mod.hist_density(g_expt, bins)
    _, ph_dens   = mod.hist_density(g_ph, bins)

    # 訓練分布が保存されていれば 95%CI と平均
    mean_dens = lo = hi = None
    try:
        import pickle, os
        if os.path.exists('gs.pkl'):
            with open('gs.pkl', 'rb') as f:
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
        ax.fill_between(c, lo, hi, color='lightgray', alpha=0.45, label='Train 95% CI')
    if mean_dens is not None:
        ax.plot(c, mean_dens, color='gray', lw=2, label='Train mean')

    ax.plot(c, ph_dens,   lw=2, label='Phantom', ls='--', c='k')
    ax.plot(c, expt_dens, lw=2, label='Experiment', ls='-', c='k')
    #ax.set_ylim(0, 0.05)
    ax.set_yscale('log')
    ax.set_xlabel(r'$|\nabla T|_{qtof95}$')
    ax.set_ylabel('Probability Density')
    #ax.legend(loc='best', frameon=False)
    ax.legend(loc='upper right', ncol=2, fontsize=8, frameon=True)
    ax.set_title(r'Probability Density of $|\nabla T|_{qtof95}$')


def draw_corr(ax):
    # （元4の図）勾配 vs π 相関
    mask    = mod.get_mask(expt=True)
    mask_ph = mod.get_mask(expt=False)
    pi_s  = (mod.get_pi().mean(axis=-1))[~mask].ravel()
    pi_ph = (mod.get_pi_phantom().mean(axis=-1))[~mask_ph].ravel()
    g_s, g_ph = mod.get_gs()

    g_min = np.nanmin([np.nanmin(g_s),  np.nanmin(g_ph)])
    g_max = np.nanmax([np.nanmax(g_s),  np.nanmax(g_ph)])
    bins = np.linspace(g_min, g_max, 41)

    xc_s, mean_s, sem_s = mod.bin_average(g_s,  pi_s,  bins)
    xc_p, mean_p, sem_p = mod.bin_average(g_ph, pi_ph, bins)

    #ax.errorbar(xc_s, mean_s, yerr=sem_s, fmt='-o', label='Experiment', capsize=3, c='tab:green')
    ax.errorbar(xc_s, mean_s, yerr=sem_s, color='black', marker='o', lw=0.7, ls='-', label='Experiment', capsize=3)
    #ax.errorbar(xc_p, mean_p, yerr=sem_p, fmt='-s', label='Phantom',    capsize=3, c='tab:orange')
    ax.errorbar(xc_p, mean_p, yerr=sem_p, color='black', marker='s', lw=0.7, ls='-', mfc='white', mew=1, label='Phantom', capsize=3)
    ax.set_xlabel(r'Transmission gradient $|\nabla T|_{qtof95}$')
    ax.set_ylabel('Uncertainty violation rate π')
    #ax.grid(True, alpha=0.3); ax.legend(loc='best', frameon=False)
    ax.set_title(r'Correlation between $|\nabla T|_{qtof95}$ and π')
    ax.legend(loc='best', frameon=False)


def draw_line(ax):
    # （元3の図）マップを描かず、回転後の行平均ラインだけを描く
    import cv2
    g2d_s, g2d_ph = mod.get_g2ds()
    H, W = g2d_s.shape
    center = (int(H/2), int(W/2))
    angle, angle_ph = 5.5, -7.0
    scale = 1.0
    trans    = cv2.getRotationMatrix2D(center, angle,    scale)
    trans_ph = cv2.getRotationMatrix2D(center, angle_ph, scale)
    g2d_s_rot  = cv2.warpAffine(g2d_s,  trans,    (W, H))
    g2d_ph_rot = cv2.warpAffine(g2d_ph, trans_ph, (W, H))[::-1, :]

    ax.plot(g2d_s_rot.mean(axis=-1), ls='-', label='Experiment', c='k')
    ax.plot(g2d_ph_rot.mean(axis=-1), ls='--', label='Phantom', c='k')
    ax.set_xlim(0, 50)  # 元コード準拠
    ax.set_xlabel('Transverse position / pixel')
    ax.set_ylabel(r'$<|\nabla T|_{qtof95}>_{longitude}$')
    ax.legend(loc='best', frameon=False)
    ax.set_title(r'Spatial Distribution of $<|\nabla T|_{qtof95}>_{longitude}$')


def get_axes_and_fig():
    fig = plt.figure()
    # ★ インチで厳密配置（必要なら値を調整）
    left = 0.65
    right = 0.45
    top = 0.35
    bottom = 0.55
    panel_w = 4.8
    h1 = 1.55
    h2 = 1.55
    h3 = 1.55
    #h4 = 1.25
    s12 = 0.28*2.5
    s23 = 0.25*2.5
    #s34 = 0.28*2
    """
    1列×4段の Axes をインチ指定で厳密配置し、[ax1, ax2, ax3, ax4] を返す。
    """
    # 図全体サイズ（inch）
    w = left + panel_w + right
    #h = top + (h1 + s12 + h2 + s23 + h3 + s34 + h4) + bottom
    h = top + (h1 + s12 + h2 + s23 + h3) + bottom
    # 図の物理サイズを決め打ち
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


def main(output="grad_all_in_one_inches.pdf"):
    fig, ax1, ax2, ax3 = get_axes_and_fig()
    # 並べ順は (1) → (2) → (4) → (3)
    draw_density(ax2)
    #draw_survival(ax3)
    draw_corr(ax1)   # ← 元「4番目」
    draw_line(ax3)   # ← 元「3番目（マップなし）」

    # 体裁：上3段の x ラベルは隠す（横位置は inch 固定なのでズレません）
    #for a in [ax1, ax2, ax3]:
    #    a.tick_params(labelbottom=False)

    #fig.savefig(output, dpi=300)
    plt.show()
    #print(f"Saved -> {output}")

if __name__ == "__main__":
    main()

