#!/usr/bin/env python
# file: demo_inch_layout_example.py
import numpy as np
import matplotlib.pyplot as plt

def get_data():
    import pickle
    with open('/home/kazu/desktop/240424/uNID_data_KO/211/openbeam.pkl',
              'rb') as f:
        openbeam = pickle.load(f)
        openbeam_noisy = pickle.load(f)
        sample = pickle.load(f)
        sample_noisy = pickle.load(f)
        etc = pickle.load(f)
    with open('/home/kazu/desktop/240424/uNID_data_KO/255/openbeamstrd.pkl',
              'rb') as f:
        openbeam_strd = pickle.load(f)
        openbeam_noisy_strd = pickle.load(f)
        sample_strd = pickle.load(f)
        sample_noisy_strd = pickle.load(f)
        etc_strd = pickle.load(f)
    transmission_noisy = sample_noisy/openbeam_noisy
    transmission = sample/openbeam*315715./553690
    transmission_noisy_strd = sample_noisy_strd/openbeam_noisy_strd
    transmission_strd = sample_strd/openbeam_strd*315715./553690
    return transmission_noisy, transmission, transmission_noisy_strd,\
        transmission_strd

def plot_expt_transmission(
    D_n: np.ndarray,
    D: np.ndarray,
    D_ns: np.ndarray,
    D_s: np.ndarray,
    cmap: str = "viridis",
) -> plt.Figure:
    #def plot_inch_layout_demo(cmap="viridis"):
    #"""
    #4列×2行のレイアウトを inch 単位で完全制御するデモ。
    #上段は画像 + カラーバー、下段は線グラフ。
    #本体軸の横幅・左位置を上下で一致させる。
    #"""
    # -------------------- 図サイズ（inch） --------------------
    fig_w_in = 12.0
    fig_h_in = 3.0
    fig = plt.figure(figsize=(fig_w_in, fig_h_in))

    # -------------------- レイアウト寸法（inch） --------------------
    ncols = 4

    # 余白・間隔（好みに合わせて調整可）
    left_in   = 0.35   # 左余白
    right_in  = 0.25   # 右余白
    bottom_in = 0.30   # 下余白
    top_in    = 0.20   # 上余白
    col_gap_in = 0.40  # 列間の隙間（本体→次本体）
    row_gap_in = 0.35  # 行間（下段→上段）

    # 本体の高さ（上段・下段）
    top_h_in    = 1.10
    bottom_h_in = 1.10

    # カラーバーの太さと本体からの隙間（inch）
    cbar_w_in   = 0.04
    cbar_pad_in = 0.02

    # -------------------- 本体の横幅（inch）を自動算出 --------------------
    main_w_in = (
        fig_w_in
        - left_in
        - right_in
        - ncols * (cbar_w_in + cbar_pad_in)
        - (ncols - 1) * col_gap_in
    ) / ncols
    if main_w_in <= 0:
        raise ValueError("寸法の合計が図幅を超えています。余白や隙間、カラーバー寸法を見直してください。")

    # inch→相対座標（0–1）変換
    W, H = fig.get_size_inches()
    xr = lambda inch: inch / W
    yr = lambda inch: inch / H

    # 上段・下段の y 位置
    y_bottom = yr(bottom_in)
    y_top    = yr(bottom_in + bottom_h_in + row_gap_in)

    # ダミーデータ（上段用 2D、下段用 1D）
    # 実データに置き換えるなら、あなたの配列に差し替えてください。
    #rng = np.random.default_rng(0)
    #img_shape = (80, 120)   # (y, x)
    #line_len  = 200

    labels = ['1/7 data ', '1/1 data', '1/7 data with binning', '1/1 data with binning']
    for didx, (data, dataname) in enumerate(zip([D_n, D, D_ns, D_s], labels)):

                                                    #for i, dataname in enumerate(labels):
        x_in = left_in + didx * (main_w_in + cbar_w_in + cbar_pad_in + col_gap_in)

        # --- 上段（画像） 本体 ---
        ax_img = fig.add_axes([xr(x_in), y_top, xr(main_w_in), yr(top_h_in)])
        #img = rng.normal(loc=0.9, scale=0.2, size=img_shape)
        #im = ax_img.imshow(data, origin="lower", cmap=cmap, vmin=0.2, vmax=1.4, aspect="auto")
        im = ax_img.imshow(
            data[100].T,
            origin="lower",
            cmap=cmap,
            aspect="auto",
            interpolation="nearest",
            vmin = 0.2,
            vmax = 1.4
        )

        ax_img.set_title(dataname, fontsize=9)
        ax_img.set_xlabel('x / ch', labelpad=1)
        ax_img.set_ylabel('y / ch', labelpad=1)
        ax_img.tick_params(length=2, labelsize=8, pad=1.1)

        # 上段（画像） カラーバー（inch固定）
        if didx == len(labels)-1:
            cax = fig.add_axes([xr(x_in + main_w_in + cbar_pad_in), y_top, xr(cbar_w_in), yr(top_h_in)])
            cbar = fig.colorbar(im, cax=cax)
            cbar.set_label("Transmission", fontsize=8, labelpad=1)
            cbar.ax.tick_params(direction='in', length=1.3, pad=1, labelsize=8)

        # --- 下段（線グラフ） 本体（上段と同じ幅＆左位置）---
        ax_line = fig.add_axes([xr(x_in), y_bottom, xr(main_w_in), yr(bottom_h_in)])
        #line = rng.normal(loc=0.9, scale=0.25, size=line_len)
        ax_line.plot(np.arange(data.shape[0])*20.+23000, data[:, 100, 36], marker='o', ms=2, lw=0)
        ax_line.set_ylim([0, 1.5])
        ax_line.set_xticks([23000, 24000, 25000, 26000])
        ax_line.set_xlabel(r'TOF / $\mu$s', labelpad=1)
        ax_line.set_ylabel('Transmission', labelpad=1)
        ax_line.tick_params(direction="in", labelsize=8, pad=2.)

    return fig


if __name__ == "__main__":
    transmission_noisy, transmission, transmission_noisy_strd,\
            transmission_strd = get_data()
    figure = plot_expt_transmission(transmission_noisy, transmission,
                                   transmission_noisy_strd, transmission_strd,
                                   cmap="viridis")
    figure.savefig("fig_expt_transmission.eps", dpi=160, bbox_inches="tight")
    plt.show()
    

