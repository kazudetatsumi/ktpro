#!/usr/bin/env python
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/kazu/ktpro")
from sf_bi_spike_at_boudary_phantom import compare_images4_2d_phantom,\
        compare_images4_d_phantom
from sf_bi_spike_at_boudary import compare_images4_2d, compare_images4_d
import matplotlib as mpl

mpl.rcParams.update({
    # EPS/PDF で Type 42（TrueType）を埋め込み → 文字がテキストとして保たれやすい
    "ps.fonttype": 42,
    "pdf.fonttype": 42,

    # Unicode/glyph が豊富で Matplotlib と相性の良いデフォルト
    "font.family": "DejaVu Sans",   # あるいは "DejaVu Serif" でもOK
    # 数式（ギリシャ文字・\AA 等）を STIX 系に揃える
    "mathtext.fontset": "stix",
    # マイナス記号を正しいUnicode（−）にする（環境で豆腐回避）
    "axes.unicode_minus": False,
})



def compare_2d():
    fig = plt.figure(figsize=(11, 10))
    subfigs = fig.subfigures(2, 1, hspace=0.0)
    for sf, plotter in zip(subfigs, [compare_images4_2d_phantom,
                                     compare_images4_2d]):
        plotter(sf)
    plt.savefig('compare_spikes_2d.pdf')
    plt.show()


def compare_1d():
    fig = plt.figure(figsize=(4, 6))
    subfigs = fig.subfigures(2, 1, hspace=0.0)
    for sf, plotter in zip(subfigs, [compare_images4_d_phantom,
                                     compare_images4_d]):
        plotter(sf)
    plt.savefig('compare_spikes_1d.eps')
    plt.show()


compare_2d()
