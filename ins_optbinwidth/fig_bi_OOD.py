#!/usr/bin/env python
# integrate_three_columns.py
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
#from plot_column import plot_one_column_in_subfigure
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
    top    = ref_top_ax.get_position().y1
    bottom = ref_bottom_ax.get_position().y0
    target_span = top - bottom
    if target_span <= 0:
        return  # 異常系は無視

    # 3) 現在の第3カラムの帯（Figure座標の 0–1）
    pos = [ax.get_position() for ax in col3_axes_top2bottom]
    cur_bottom = min(p.y0 for p in pos)
    cur_top    = max(p.y1 for p in pos)
    cur_span   = cur_top - cur_bottom
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
PANEL_W_IN  = 3.0
MAP_H_IN    = 2.2
LINE_H_IN   = 1.2
LM, RM, TM, BM = 0.4, 0.4, 0.25, 0.45  # 外側（図全体）の余白（インチ）
WSPACE_IN   = 0.5                      # 列間の見かけ隙間（インチ相当）
HSPACE_IN   = 0.3                      # 各列内部の行間は列関数側で調整

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

# ====== データの読み込み（例） ======
# 実データに置換してください
#def load_data_A():
#    # 形・意味はユーザーの既存配列に揃える
#    cond = np.random.rand(64, 64) + 1e-3
#    maxpos = (30, 40)
#    maxposin = (28, 35)
#    poslists = [(20, 22)]
#    maxstd = np.random.rand(64, 64, 100)
#    std_phantom_sample = np.random.rand(64, 64, 100)
#    return cond, maxpos, maxposin, poslists, maxstd, std_phantom_sample
#
#dataA = load_data_A()
#dataB = load_data_A()
#dataC = load_data_A()
#
# ====== 各列に描画 ======
#for sf, data, label in zip(subfigs, [dataA, dataB, dataC], ["(A)", "(B)", "(C)"]):
for sf, plotter in zip(subfigs, [check_std, check_std_phantom, draw_grad_column_in_sf]):
    axes = plotter(
        sf,
    )
    # 列タイトル（subfigure 内の最上段に）
    # axes["ax_map"].set_title(f"{label} Column title", loc="left", pad=6)



# --- 3カラムを描画した直後に追記 ---
# 基準はカラムA（subfigs[0]）を採用（BでもOK）
ref_top_ax    = find_top_map_ax_in_sf(subfigs[0])         # A列のマップ上端
ref_bottom_ax = find_bottom_plot_ax_in_sf(subfigs[0])     # A列の最下段プロット下端

# 第3カラムの Axes（上→下）
#col3 = draw_grad_column_in_sf(subfigs[2], title_left="(C)")
ax_top, ax2, ax3, ax_bottom = axes

# 縦帯（top/bottom）をインチ基準で一致させる
align_col3_to_band_in_inches(fig, subfigs[2], [ax_top, ax2, ax3, ax_bottom],
                             ref_top_ax, ref_bottom_ax)

plt.show()

# ====== 仕上げ ======
#mpl.rcParams["pdf.fonttype"] = 42
#mpl.rcParams["ps.fonttype"] = 42
#fig.savefig("three_columns.pdf", dpi=300)
# fig.savefig("three_columns.eps", dpi=300)  # EPSが必要ならこちら

