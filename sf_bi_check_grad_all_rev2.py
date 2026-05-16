#!/usr/bin/env/python
# integrate_grad_onefig.py
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec



# --- (新規) パネル(1)〜(3)で共通に使う xlim 計算 ---
def get_common_xlim():
    # ΔT の最小・最大を実験/ファントムから取得
    g_expt, g_ph = mod.get_gs()
    xmin = np.nanmin([np.nanmin(g_expt), np.nanmin(g_ph)])
    xmax = np.nanmax([np.nanmax(g_expt), np.nanmax(g_ph)])
    if not np.isfinite(xmin) or not np.isfinite(xmax):
        return None
    pad = 0.02 * (xmax - xmin) if xmax > xmin else 1.0  # 2%だけ余白
    return (xmin - pad, xmax + pad)


def draw_grad_column_in_sf(sf,
                           fig_width_in=6.6,
                           h_each_in=1.6,           # (1)〜(3) 各段の高さ
                           gap_between_blocks_in=1.6,  # (3) と (4) の間隔（ここを調整すれば「少し間」）
                           h4_in=1.2,                 # (4) の高さ
                           top_margin_in=0.3,
                           bottom_margin_in=0.4,
                           left_margin_in=0.6,
                           right_margin_in=0.6,
                           ):
    mpl.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42, 
        "font.size": 9, "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "legend.fontsize": 8,
    })  

    fig_height_in = (top_margin_in + 3*h_each_in
                     + gap_between_blocks_in + h4_in + bottom_margin_in)
    #fig = plt.figure(figsize=(fig_width_in, fig_height_in))

    # ===== 2段グリッド：上段＝(1)-(3)ブロック, 下段＝(4) =====
    gs = sf.add_gridspec(
        2, 1,
        height_ratios=[3*h_each_in, h4_in],
        left=left_margin_in/fig_width_in,
        right=1 - right_margin_in/fig_width_in,
        bottom=bottom_margin_in/fig_height_in,
        top=1 - top_margin_in/fig_height_in,
        # (3) と (4) の間にだけ隙間を入れる
        #hspace=gap_between_blocks_in/fig_height_in
        hspace=0.4
    )

    # ===== 上段ブロック： (1)〜(3) を “横軸共有・上下隙間ゼロ” で縦積み =====
    area = sf.add_subplot(gs[0, 0])
    area_pos = area.get_position()   # 図座標
    sf.delaxes(area)                # ダミーなので削除

    h_each = area_pos.height / 3.0
    left = area_pos.x0
    width = area_pos.width
    bottom3 = area_pos.y0
    bottom2 = bottom3 + h_each
    bottom1 = bottom2 + h_each

    # まず最下段 ax3 を作って、それを sharex に指定
    ax3 = sf.add_axes([left, bottom3, width, h_each])           # (3) Correlation
    ax2 = sf.add_axes([left, bottom2, width, h_each], sharex=ax3)  # (2) Survival
    ax1 = sf.add_axes([left, bottom1, width, h_each], sharex=ax3)  # (1) Hist

    # 描画（既存関数をそのまま使用）
    draw_panel1_hist(ax1)
    draw_panel2_survival(ax2)
    draw_panel3_corr(ax3)

    # (1)〜(3) の xlim/xticks を共通化
    xlim = get_common_xlim()
    if xlim is not None:
        for ax in (ax1, ax2, ax3):
            ax.set_xlim(xlim)

    # 上2段は x ラベル・x 目盛りを消して “くっつけた” 見た目に
    for ax in (ax1, ax2):
        ax.set_xlabel('')
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.tick_params(labelbottom=False)

    # ===== 下段ブロック： (4) は横軸が別なので独立の軸 =====
    ax4 = sf.add_subplot(gs[1, 0])
    draw_panel4_line(ax4)

    # 左右端の完全一致
    align_left_right(sf, [ax1, ax2, ax3, ax4])
    for ax in (ax1, ax2, ax3, ax4):
        ax.tick_params(direction='in', which='both', top=True, right=True)

    # 保存と表示（元のファイル名を踏襲）
    #fig.savefig("grad_all_in_one.pdf", dpi=300)
    #plt.show()
    return [ax1, ax2, ax3, ax4]




# 同ディレクトリの既存スクリプトをモジュールとして使う
import fig_bi_check_grad as mod   # ← あなたの元ファイル名 [1](https://jaeagojp-my.sharepoint.com/personal/tatsumi_kazuyoshi_jaea_go_jp/Documents/Microsoft%20Copilot%20Chat%20%E3%83%95%E3%82%A1%E3%82%A4%E3%83%AB/fig_bi_check_grad.py)

# ---- 1) パネル1：勾配ヒストグラム（Train帯はgs.pklがあれば描く） ----
def draw_panel1_hist(ax):
    # 既存ヘルパを流用
    g_expt, g_ph = mod.get_gs()  # 実験/模擬刀の ∇T flatten 配列を返す前処理 [1](https://jaeagojp-my.sharepoint.com/personal/tatsumi_kazuyoshi_jaea_go_jp/Documents/Microsoft%20Copilot%20Chat%20%E3%83%95%E3%82%A1%E3%82%A4%E3%83%AB/fig_bi_check_grad.py)
    all_vals = np.concatenate([g_expt, g_ph])
    eps = 1e-8
    vmin = np.percentile(all_vals, 0.1)
    vmax = np.percentile(all_vals, 99.9)
    bins = np.exp(np.linspace(np.log(max(vmin, eps)), np.log(vmax), 30))
    c, expt_dens = mod.hist_density(g_expt, bins)
    _, ph_dens   = mod.hist_density(g_ph, bins)
    # 可能なら訓練分布の帯（gs.pklの3本目）を描く
    mean_dens = lo = hi = None
    try:
        # compare_grads() と同じファイル名に合わせて読む [1](https://jaeagojp-my.sharepoint.com/personal/tatsumi_kazuyoshi_jaea_go_jp/Documents/Microsoft%20Copilot%20Chat%20%E3%83%95%E3%82%A1%E3%82%A4%E3%83%AB/fig_bi_check_grad.py)
        import pickle, os
        if os.path.exists('gs.pkl'):
            with open('gs.pkl', 'rb') as f:
                _ = pickle.load(f)            # g_sample_expt
                _ = pickle.load(f)            # g_sample_phantom
                g_x_train_noisy = pickle.load(f)
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
        ax.fill_between(c, lo, hi, color='tab:blue', alpha=0.15, label='Train 95% CI')
    if mean_dens is not None:
        ax.plot(c, mean_dens, color='tab:blue', lw=2, label='Train mean')

    ax.plot(c, ph_dens,   color='tab:orange', lw=2, label='Phantom')
    ax.plot(c, expt_dens, color='tab:green',  lw=2, label='Experiment')
    ax.set_ylim(0, 0.05)
    ax.set_xlabel('∇T')
    ax.set_ylabel('Density')
    ax.legend(loc='best', frameon=False)
    #ax.set_title('(1) Gradient histogram')

# ---- 2) パネル2：サバイバル曲線（gs.pkl に train が無い場合は Expt/Ph のみ） ----
def draw_panel2_survival(ax):
    # expt / phantom から上側確率S(g)=P(G>g)を描く
    g_expt, g_ph = mod.get_gs()  # [1](https://jaeagojp-my.sharepoint.com/personal/tatsumi_kazuyoshi_jaea_go_jp/Documents/Microsoft%20Copilot%20Chat%20%E3%83%95%E3%82%A1%E3%82%A4%E3%83%AB/fig_bi_check_grad.py)
    # 可能なら train も
    g_train_all = None
    try:
        import pickle, os
        if os.path.exists('gs.pkl'):
            with open('gs.pkl', 'rb') as f:
                _ = pickle.load(f); _ = pickle.load(f); g_x_train_noisy = pickle.load(f)
            g_train_all = np.concatenate(g_x_train_noisy)
    except Exception:
        pass

    # グリッドは元コードの q95 近辺に倣って適当な範囲を組む
    def surv(x, grid=None):
        x = np.sort(x[np.isfinite(x)])
        if grid is None:
            grid = np.quantile(x, np.linspace(0, 0.999, 200))
        S = 1.0 - np.searchsorted(x, grid, side='right')/x.size
        return grid, S

    if g_train_all is not None:
        q95 = np.quantile(g_train_all, 0.95)
        q99 = np.quantile(g_train_all, 0.99)
        #grid = np.linspace(q95, max(np.nanmax(g_expt), np.nanmax(g_ph), np.nanmax(g_train_all)), 200)
        grid = np.linspace(min(np.nanmin(g_expt), np.nanmin(g_ph), np.nanmin(g_train_all)), max(np.nanmax(g_expt), np.nanmax(g_ph), np.nanmax(g_train_all)), 200)
        xt, St = surv(g_train_all, grid)
        ax.semilogy(xt, St, color='tab:blue', label='Train (all)')
        ax.axvline(q95, color='k', ls='--', alpha=0.5); ax.axvline(q99, color='k', ls=':', alpha=0.5)
    else:
        grid = None

    xm, Sm = surv(g_ph,  grid)
    xe, Se = surv(g_expt, grid)
    ax.semilogy(xm, Sm, color='tab:orange', label='Phantom')
    ax.semilogy(xe, Se, color='tab:green',  label='Experiment')
    ax.set_xlabel('∇T'); ax.set_ylabel("Survival S(g)=P(∇T>g)")
    tmp_ylim = (ax.get_ylim())
    ax.set_ylim([tmp_ylim[0], tmp_ylim[1]*2])
    ax.legend(loc='best', frameon=False)
    #ax.set_title('(2) Survival curves')

# ---- 3) パネル3：「勾配 vs π」相関（4番目の図） ----
def draw_panel3_corr(ax):
    # compare_grads_and_pi() 末尾の相関エラーバー部を移植 [1](https://jaeagojp-my.sharepoint.com/personal/tatsumi_kazuyoshi_jaea_go_jp/Documents/Microsoft%20Copilot%20Chat%20%E3%83%95%E3%82%A1%E3%82%A4%E3%83%AB/fig_bi_check_grad.py)
    mask = mod.get_mask(expt=True)
    mask_ph = mod.get_mask(expt=False)
    # π（過小評価発生率）と ∇T を flatten
    pi_s  = (mod.get_pi().mean(axis=-1))[~mask].ravel()
    pi_ph = (mod.get_pi_phantom().mean(axis=-1))[~mask_ph].ravel()
    g_s, g_ph = mod.get_gs()
    # 共通ビン
    g_min = np.nanmin([np.nanmin(g_s),  np.nanmin(g_ph)])
    g_max = np.nanmax([np.nanmax(g_s),  np.nanmax(g_ph)])
    bins = np.linspace(g_min, g_max, 41)
    # ビン平均と標準誤差
    xc_s, mean_s, sem_s = mod.bin_average(g_s,  pi_s,  bins)
    xc_p, mean_p, sem_p = mod.bin_average(g_ph, pi_ph, bins)

    ax.errorbar(xc_s, mean_s, yerr=sem_s, fmt='-o', label='Experiment', capsize=3)
    ax.errorbar(xc_p, mean_p, yerr=sem_p, fmt='-s', label='Phantom',    capsize=3)
    ax.set_xlabel('$Transmission\ gradient\ ∇T\ /\ ch^{-1}$')
    ax.set_ylabel('Underestimation rate π')
    ax.grid(True, alpha=0.3); ax.legend(loc='best', frameon=False)
    tmp_ylim = (ax.get_ylim())
    ax.set_ylim([tmp_ylim[0], tmp_ylim[1]+0.1])
    #ax.set_title('(3) Correlation: gradient vs π')

# ---- 4) パネル4：マップ無し・下段ラインだけ（3番目の“下段”を抽出） ----
def draw_panel4_line(ax):
    # compare_grads_and_pi() 中の「下段ライン」だけを再現（マップは描かない） [1](https://jaeagojp-my.sharepoint.com/personal/tatsumi_kazuyoshi_jaea_go_jp/Documents/Microsoft%20Copilot%20Chat%20%E3%83%95%E3%82%A1%E3%82%A4%E3%83%AB/fig_bi_check_grad.py)
    import cv2
    g2d_s, g2d_ph = mod.get_g2ds()  # 2Dの ∇T マップ（expt/phantom） [1](https://jaeagojp-my.sharepoint.com/personal/tatsumi_kazuyoshi_jaea_go_jp/Documents/Microsoft%20Copilot%20Chat%20%E3%83%95%E3%82%A1%E3%82%A4%E3%83%AB/fig_bi_check_grad.py)
    H, W = g2d_s.shape
    center = (int(H/2), int(W/2))
    angle, angle_ph = 5.5, -7.0
    scale = 1.0
    trans     = cv2.getRotationMatrix2D(center, angle,    scale)
    trans_ph  = cv2.getRotationMatrix2D(center, angle_ph, scale)
    g2d_s_rot  = cv2.warpAffine(g2d_s,  trans,    (W, H))
    g2d_ph_rot = cv2.warpAffine(g2d_ph, trans_ph, (W, H))[::-1, :]

    ax.plot(g2d_s_rot.mean(axis=-1),  label='Experiment')
    ax.plot(g2d_ph_rot.mean(axis=-1), label='Phantom')
    ax.set_xlim(0, 50)  # 元コード準拠
    ax.set_xlabel('y\' / ch')
    ax.set_ylabel('Mean ∇T over x\'')
    ax.legend(loc='best', frameon=False)
    #ax.set_title('(4) Ridge-aligned line profile (maps omitted)')

# ---- 左右端を“機械的に”完全一致させる（post-draw整列） ----
def align_left_right(fig, axes):
    # constrained_layout などの再配置が終わった後で呼ぶ
    fig.canvas.draw()
    for a in axes: a.set_in_layout(True)
    pos = [a.get_position() for a in axes]
    left  = max(p.x0 for p in pos)   # もっとも右に寄っている左端
    right = min(p.x1 for p in pos)   # もっとも左に寄っている右端
    for a, p in zip(axes, pos):
        a.set_position([left, p.y0, right - left, p.height])


#if __name__ == "__main__":
#    main()
