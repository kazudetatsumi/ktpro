#!/usr/bin/env python
# This script forms a figure containing three columns: The first column
# shows distributions of the std of estimated neutron counts in the exp data
# point, the second those in the phantom data point.The thrid column frst row
# compares histograms of transmission gradients in the experimental, phantom
# and training data sets.
# The third column second row compares survival rates of transmission gradients
# in the experimental, phantom and training data sets.
# The third column third row shows the dependence of the probablitiy of the
# std of the estimated neutron coutns in the experimental or phantom data
# points exceeds the maximum std of another set of training data points on the
# transmission gradients. The third column fourth row shows the transmission
# gradients averaged over the y' axis are plotted with respect to the x' axis,
# where the (x, y) was rotated to (x', y') so that the sword tou is parallel to
# the x'.
# To prepare a paper manuscript on denoising the Bragg Edge image  data.
# Kazuyoshi TATSUMI 2026/01/23
from pathlib import Path
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/kazu/ktpro")
from sf_bi_check_std import check_std
from sf_bi_check_std_phantom import check_std_phantom
from sf_bi_check_grad_all_rev2 import draw_grad_column_in_sf


# --- ユーティリティ：SubFigure内の「最上段マップ」と「最下段ライン」を探す（返り値: Axes） ---
def find_top_map_ax_in_sf(sf):
    # imshow 等を持つ Axes のうち、最も上にあるもの
    cands = [ax for ax in sf.axes if getattr(ax, "images", None)]
    if not cands:
        # 画像が無い場合のフォールバック：単に最上段の Axes
        cands = list(sf.axes)
    return max(cands, key=lambda ax: ax.get_position().y1)


def find_bottom_plot_ax_in_sf(sf):
    # 折れ線等を持つ Axes のうち、最も下にあるもの
    cands = [ax for ax in sf.axes if getattr(ax, "lines", None)]
    if not cands:
        cands = list(sf.axes)
    return min(cands, key=lambda ax: ax.get_position().y0)


def align_col3_to_band_in_inches(fig, sf_col3, col3_axes_top2bottom,
                                 ref_top_ax, ref_bottom_ax):
    """
    fig                : Figure
    sf_col3            : 第3カラムの SubFigure（未使用だがシグネチャは互換のまま）
    col3_axes_top2bottom : 第3カラムの Axes [ax_top, ax2, ax3, ax_bottom]（上→下）
    ref_top_ax         : 参照カラムの「マップ軸」（上端 y1 を合わせる）
    ref_bottom_ax      : 参照カラムの「最下段プロット軸」（下端 y0 を合わせる）
    """
    # 1) 最終レイアウト確定（constrained_layout=True のため必須）
    fig.canvas.draw()

    # 2) 参照の帯（Figure座標の 0–1）
    top = ref_top_ax.get_position().y1
    bottom = ref_bottom_ax.get_position().y0
    target_span = top - bottom
    if target_span <= 0:
        return  # 異常系は無視

    # 3) 現在の第3カラムの帯（Figure座標の 0–1）
    pos = [ax.get_position() for ax in col3_axes_top2bottom]
    cur_bottom = min(p.y0 for p in pos)
    cur_top = max(p.y1 for p in pos)
    cur_span = cur_top - cur_bottom
    if cur_span <= 0:
        return

    # 4) 各軸の y0/y1 を線形写像で [cur_bottom,cur_top] → [bottom,top] に再配置（xはそのまま）
    for ax, p in zip(col3_axes_top2bottom, pos):
        y0n = (p.y0 - cur_bottom) / cur_span
        y1n = (p.y1 - cur_bottom) / cur_span
        new_y0 = bottom + y0n * target_span
        new_y1 = bottom + y1n * target_span
        ax.set_position([p.x0, new_y0, p.width, new_y1 - new_y0])


# ====== 全体寸法（inch） ======
NCOLS = 3
# 各カラムの想定パネル寸法は “列内部の関数”で決めているが、
# 外側の見かけギャップや外余白を合わせるため、ここで全体サイズを設計
PANEL_W_IN = 3.0
MAP_H_IN = 2.2
LINE_H_IN = 1.2
LM, RM, TM, BM = 0.4, 0.4, 0.25, 0.45  # 外側（図全体）の余白（インチ）
WSPACE_IN = 0.5                      # 列間の見かけ隙間（インチ相当）
HSPACE_IN = 0.3                      # 各列内部の行間は列関数側で調整

# 図全体の高さは列内部と整合するように（例：上の寸法に合わせて）
fig_w = LM + NCOLS*PANEL_W_IN + (NCOLS-1)*WSPACE_IN + RM
fig_h = TM + (MAP_H_IN + 3*LINE_H_IN) + BM

# ====== Figure 作成（subfigures） ======
fig = plt.figure(figsize=(fig_w, fig_h), constrained_layout=True)
# constrained_layout のパディングをインチ感覚で設定
fig.set_constrained_layout_pads(
    w_pad=LM, h_pad=BM,
    wspace=WSPACE_IN / PANEL_W_IN,   # 列間をパネル幅比で調整
    hspace=HSPACE_IN / LINE_H_IN     # （全体の縦詰め微調整）
)

subfigs = fig.subfigures(1, NCOLS, wspace=0.02)

for sf, plotter in zip(subfigs, [check_std_phantom, check_std,
                                 draw_grad_column_in_sf]):
    axes = plotter(sf)
ref_top_ax = find_top_map_ax_in_sf(subfigs[0])
ref_bottom_ax = find_bottom_plot_ax_in_sf(subfigs[0])
align_col3_to_band_in_inches(fig, subfigs[2], axes, ref_top_ax, ref_bottom_ax)
#plt.savefig('check_std_all.eps')
plt.show()
