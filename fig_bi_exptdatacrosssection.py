#!/usr/bin/env python
# file: demo_inch_layout_example.py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 9/2*1.3, "axes.labelsize": 11/2*1.3, "xtick.labelsize": 9/2*1.3, "ytick.labelsize": 9/2*1.3,
    "legend.fontsize": 9/2*1.3,
})

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
        transmission_strd, sample_noisy, sample, sample_noisy_strd,\
        sample_strd, etc, etc_strd


def plot_expt_transmission(
    D_n: np.ndarray,
    D: np.ndarray,
    D_ns: np.ndarray,
    D_s: np.ndarray,
    dtype: str = 'Transmission',
    cmap: str = "gray_r",
) -> plt.Figure:
    fig_w = 12.0*0.65
    fig_h = 3.0*0.65
    fig = plt.figure(figsize=(fig_w, fig_h))
    ncols = 4
    left = 0.35*0.65  # left margin
    right = 0.25*0.65  # right margin
    bottom = 0.30*0.65   # bottom margin
    col_gap = 0.85*0.65  # column space
    row_gap = 0.35*0.65  # row space
    # height of top and bottom rows
    top_h = 1.10*0.65
    bottom_h = 1.10*0.65
    # color bar width and pad
    cbar_w = 0.04*0.7
    cbar_pad = 0.02*0.8

    # -------------------- width of main body -------------------
    main_w = (
        fig_w
        - left
        - right
        - ncols * (cbar_w + cbar_pad)
        - (ncols - 1) * col_gap
    ) / ncols
    if main_w <= 0:
        raise ValueError("total dimension exceeds figure width." +
                         "Reconsider your margin, color bar size")

    # inch-> fractional coordination
    W, H = fig.get_size_inches()
    # y position of top and bottom rows
    y_bottom = bottom/H
    y_top = (bottom + bottom_h + row_gap)/H
    labels = ['1/7 data ', '1/1 data', '1/7 data with binning',
              '1/1 data with binning']
    px = 100
    py = 36
    for didx, (data, dataname) in enumerate(zip([D_n, D, D_ns, D_s], labels)):
        x = left + didx * (main_w + cbar_w + cbar_pad + col_gap)
        # --- top row of images ---
        ax_img = fig.add_axes([x/W, y_top, main_w/W, top_h/H])
        im = ax_img.imshow(
            data[100].T,
            origin="lower",
            cmap=cmap,
            aspect="auto",
            interpolation="nearest",
            #vmin=0.2,
            #vmax=1.4,
        )
        ax_img.set_title(dataname, fontsize=9*0.65)
        ax_img.set_xlabel('x / ch', labelpad=1*0.65)
        ax_img.set_ylabel('y / ch', labelpad=1*0.65)
        ax_img.annotate('', xy=(px, py), xytext=(px+9, py+9),
                        textcoords='data', color='white',
                        arrowprops=dict(arrowstyle="->", color='k'))
        if didx >= 2:
            ax_img.tick_params(direction='in', labelsize=8*0.65, pad=2,
                               color='w')
        else:
            ax_img.tick_params(direction='in', labelsize=8*0.65, pad=2,
                               color='k')
        # top row color bars
        if ((dtype == 'Transmission') and (didx == ncols-1)) \
                or dtype == 'Neutron count':
            cax = fig.add_axes([(x + main_w + cbar_pad)/W, y_top, cbar_w/W,
                                top_h/H])
            cbar = fig.colorbar(im, cax=cax)
            cbar.set_label(dtype, labelpad=1*0.65)
            cbar.ax.tick_params(direction='in', length=1.3*0.8, pad=1*0.8,)
        # --- bottom row ---
        ax_line = fig.add_axes([x/W, y_bottom, main_w/W, bottom_h/H])
        ax_line.plot(np.arange(data.shape[0])*20.+23000,
                     data[:, px, py], marker='o', ms=2*0.65, lw=0, c='k')
        #ax_line.set_ylim([0, 1.5])
        ax_line.set_xticks([23000, 24000, 25000, 26000])
        ax_line.set_xlabel(r'TOF / $\mu$s', labelpad=2)
        ax_line.set_ylabel(dtype, labelpad=2*0.65)
        ax_line.tick_params(direction="in", labelsize=8*0.65, pad=2)
    return fig


if __name__ == "__main__":
    transmission_noisy, transmission, transmission_noisy_strd,\
        transmission_strd, sample_noisy, sample, sample_noisy_strd,\
        sample_strd, etc, etc_strd = get_data()
    """
    figure = plot_expt_transmission(transmission_noisy, transmission,
                                    transmission_noisy_strd, transmission_strd,
                                    dtype='Transmission',
                                    cmap="gray")
    """
    figure = plot_expt_transmission(sample_noisy, sample,
                                    sample_noisy_strd, sample_strd,
                                    dtype='Neutron count',
                                    cmap="gray_r")
    figure.savefig("fig_expt_neutroncount.eps", dpi=160, bbox_inches="tight")
    plt.show()
