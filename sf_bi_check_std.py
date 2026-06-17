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




def get_tdata():
    tdatafile = 'tdata.pkl'
    if not os.path.exists(tdatafile):
        for i in range(0, 5):
            with open('./seed' + str(i) + '/valtesttot.pkl', 'rb') as f:
                data = pickle.load(f)
                var = pickle.load(f)
            if i == 0:
                tdata = data.squeeze()[np.newaxis]
                tvar = var.squeeze()[np.newaxis]
            else:
                tdata = np.vstack((tdata, data.squeeze()[np.newaxis]))
                tvar = np.vstack((tvar, var.squeeze()[np.newaxis]))
        with open('tdata.pkl', 'wb') as f:
            pickle.dump(tdata.astype('float32'), f, 4)
            pickle.dump(tvar.astype('float32'), f, 4)
    else:
        print('reading ' + tdatafile)
        with open(tdatafile, 'rb') as f:
            tdata = pickle.load(f)
            tvar = pickle.load(f)
    return tdata, tvar


def get_mask(transmission):
    import cv2
    return cv2.GaussianBlur(transmission.sum(axis=-1), (5, 5), 0) > 236


def check_std(sf,
    # ===== レイアウト制御（インチ） =====
    axes_width_in=3.0,      # 1列目（マップ/各グラフ）の軸の横幅 [inch]
    map_height_in=2.2,      # 1行目（マップ）の縦幅 [inch]
    line_height_in=1.2,     # 2〜4行目の各グラフの縦幅 [inch]
    left_margin_in=0.6,     # 左余白 [inch]
    right_margin_in=0.6,    # 右余白 [inch]（cbar より右の余白）
    top_margin_in=0.3,      # 上余白 [inch]
    bottom_margin_in=0.4,   # 下余白 [inch]
    # ===== カラーバー（インチ） =====
    cbar_width_in=0.05,     # カラーバー幅 [inch]
    cbar_pad_in=0.02,       # マップとカラーバーの水平距離 [inch]
    # ===== 行間の空き（マップ行とグラフブロックの間、figure 座標で指定） =====
    hspace=0.0
):
    # ---- データ前処理（既存ロジック） ----
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

    _mask = get_mask(transmission)
    mask = get_mask(transmission)
    for iy in range(_mask.shape[1]):
        for ix in range(_mask.shape[0]):
            if (mask[ix-3, iy]) & (not mask[ix, iy]):
                _mask[ix, iy] = True

    poslists = [[40, 130], [30, 50], [30, 30], [40, 15], [66, 79]]
    cond = (maxstd < std_sample).sum(axis=-1) + 1
    condin = (maxstd < std_sample).sum(axis=-1)
    cond[mask] = 2
    condin[_mask] = 0.
    maxpos = np.unravel_index(np.argmax(cond, axis=None), cond.shape)
    maxposin = np.unravel_index(np.argmax(condin, axis=None), cond.shape)

    # ===============================
    #   図の大きさ（インチ）
    #   ※ カラーバー分の横幅も図サイズに含める
    # ===============================
    fig_height_in = top_margin_in + map_height_in + 3*line_height_in + bottom_margin_in
    fig_width_in  = left_margin_in + axes_width_in + cbar_pad_in + cbar_width_in + right_margin_in

    #fig = plt.figure(figsize=(fig_width_in, fig_height_in))

    # ===============================
    #   GridSpec（2行×1列）
    #   行1: マップ, 行2: グラフ領域（後で3軸に分割）
    #   右端はカラーバー用に "cbar_pad + cbar_width" を空ける
    # ===============================
    height_ratios = [map_height_in, 3*line_height_in]
    gs = sf.add_gridspec(
        2, 1, height_ratios=height_ratios,
        left=left_margin_in/fig_width_in,
        right=1 - (right_margin_in + cbar_pad_in + cbar_width_in)/fig_width_in,
        bottom=bottom_margin_in/fig_height_in,
        top=1 - top_margin_in/fig_height_in,
        hspace=hspace
    )

    # ===== 1段目：マップ =====
    ax_map = sf.add_subplot(gs[0, 0])
    cim = ax_map.imshow(cond, norm=colors.LogNorm(vmin=cond.min(), vmax=cond.max()))

    # アノテーション
    for pidx, pos in enumerate([maxpos, maxposin, poslists[0]]):
        ax_map.annotate('#'+str(pidx+1), xy=(pos[1], pos[0]),
                        xytext=(pos[1]+5, pos[0]+5), textcoords='data',
                        color='white', fontsize=10,
                        arrowprops=dict(arrowstyle="->", color='white'))

    ax_map.set_ylabel('y / ch')
    ax_map.set_xlabel('x / ch')
    # tick を in、上下右も表示
    ax_map.tick_params(length=2, labelsize=8, pad=1.1)

    # ===== マップの右に、高さ一致の cbar を配置（add_axes で厳密配置） =====
    #map_pos = ax_map.get_position()  # figure 座標
    #cbar_width_fig = cbar_width_in / fig_width_in
    #cbar_pad_fig   = cbar_pad_in   / fig_width_in
    #cax = sf.add_axes([
    #    map_pos.x1 + cbar_pad_fig,  # left
    #    map_pos.y0,                 # bottom（マップと一致）
    #    cbar_width_fig,             # width
    #    map_pos.height              # height（マップと完全一致）
    #])
    #cbar = sf.colorbar(cim, cax=cax)

    from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
    divider = make_axes_locatable(ax_map)
    cax = divider.append_axes(
        "right",
        size=axes_size.Fixed(cbar_width_in),
        pad=axes_size.Fixed(cbar_pad_in)
    )
    cbar = sf.colorbar(cim, cax=cax)
    #ax_map.set_aspect('auto')

    cbar.set_label("# of violated\nvoxels + 1")
    cax.tick_params(direction='in', labelsize=8)

    # ===== 2段目：グラフ領域（3つの個別軸を上下隙間ゼロで配置） =====
    ax_area = sf.add_subplot(gs[1, 0])
    area_pos = ax_area.get_position()  # figure 座標での矩形
    sf.delaxes(ax_area)               # 使わないので削除

    # 3等分して上下に隙間ゼロで配置
    h_each = area_pos.height / 3.0
    left = area_pos.x0
    width = area_pos.width
    bottom4 = area_pos.y0
    bottom3 = bottom4 + h_each
    bottom2 = bottom3 + h_each

    # 最下段を基準に share するため、先に ax4 を作成
    ax4 = sf.add_axes([left, bottom4, width, h_each])
    ax3 = sf.add_axes([left, bottom3, width, h_each], sharex=ax4)
    ax2 = sf.add_axes([left, bottom2, width, h_each], sharex=ax4)

    # それぞれのプロット（従来の 2〜4 行目に対応）
    x = np.arange(maxstd.shape[-1])*20 + 23000.
    for axi, pos in zip([ax2, ax3, ax4], [maxpos, maxposin, poslists[0]]):
        axi.plot(x, maxstd[pos[0], pos[1]], label='maxstd')
        axi.plot(x, std_sample[pos[0], pos[1]], label='std_sample')
        axi.tick_params(direction='in', top=True, right=True)
        axi.set_ylim([0., 7])

    # 凡例（必要に応じて変更）
    ax2.legend(loc='best', fontsize=9)
    ax3.legend(loc='best', fontsize=9)
    ax4.legend(loc='best', fontsize=9)

    # 上2つの x 軸はラベル・目盛りを非表示（連結見せ）
    for axi in [ax2, ax3]:
        plt.setp(axi.get_xticklabels(), visible=False)
        axi.tick_params(labelbottom=False)
        axi.tick_params(labelsize=8)

    # 最下段のみ x ラベルとタイトル
    ax4.set_xlabel(r'TOF / $\mu$s')
    ax4.tick_params(labelsize=8)
    #ax4.set_title('STD comparison (rows 2–4 visually concatenated)')
    ax2.set_ylabel('std')
    ax3.set_ylabel('std')
    ax4.set_ylabel('std')
    #plt.savefig('check_std.png')
    #plt.savefig('check_std.eps')

    #plt.show()


    # --- 左右端の強制一致（post-draw） ---
    fig = sf.figure
    fig.canvas.draw()
    axes_main = [ax_map, ax2, ax3, ax4]
    for ax in axes_main:
        ax.set_in_layout(True)
    cax.set_in_layout(False)
    pos = [ax.get_position() for ax in axes_main]
    left  = max(p.x0 for p in pos)
    right = min(p.x1 for p in pos)
    for ax, p in zip(axes_main, pos):
        ax.set_position([left, p.y0, right - left, p.height])

def check_exptspectra():
    tdata = get_tdata()
    print(tdata.shape)
    for i in range(20):
        icol = np.random.randint(tdata.shape[2])
        irow = np.random.randint(tdata.shape[3])
        plt.plot(tdata[0, -2, icol, irow]/tdata[0, -1, icol, irow]*315715./553690)
        plt.savefig('exptspectra_samples.png')



