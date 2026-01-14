import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl

# Make Arial the default font for all text
mpl.rcParams['font.family'] = 'Arial'         # primary family
mpl.rcParams['font.sans-serif'] = ['Arial']   # optional: ensure sans-serif points to Arial

# Optional: math text style (see section 4)
mpl.rcParams['mathtext.fontset'] = 'dejavusans'  # keeps math readable with a sans style


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
    """Plot three aligned 2D maps (D-sum, P1, P2) with colorbars.

    Layout (no 1D profiles):
      - Top   : 2D intensity map of D (sum over z) + colorbar
      - Middle: 2D map of P1 + colorbar (shares axes with D)
      - Bottom: 2D map of P2 + colorbar (shares axes with D)

    Parameters
    ----------
    D : np.ndarray
        3D tensor with shape (Nx, Ny, Nz), axes (x, y, z).
    P1 : np.ndarray
        2D tensor with shape (Nx, Ny), axes (x, y).
    P2 : np.ndarray
        2D tensor with shape (Nx, Ny), axes (x, y).
    cmap : str, optional
        Matplotlib colormap name, by default "viridis".

    Returns
    -------
    matplotlib.figure.Figure
        The rendered figure.
    """
    print(D_n.shape, D.shape, D_ns.shape, D_s.shape)
    # Figure and grid (3 rows, 1 column)
    fig = plt.figure(figsize=(24, 6))
    grid = fig.add_gridspec(
        nrows=2,
        ncols=4,
        height_ratios=[1.0, 1.5],
        width_ratios=[1.0, 1.0, 1.0, 1.0],
        hspace=0.35,
        wspace=0.75,
    )

    for didx, (data, dataname) in enumerate(zip([D_n, D, D_ns, D_s],
                                                ['1/7 data ', '1/1 data', '1/7 data with binning',
                                                  '1/1 data with binning'])):
    # Top: D map + colorbar
        ax = fig.add_subplot(grid[0, didx])
        im = ax.imshow(
            data[100].T,
            origin="lower",
            cmap=cmap,
            aspect="auto",
            interpolation="nearest",
            vmin = 0.2,
            vmax = 1.4
        )
        ax.set_ylabel('y / ch', fontsize=18)
        ax.set_xlabel('x / ch', fontsize=18)
        ax.set_title(dataname, fontsize=18)
        cbar_d = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar_d.set_label("Transmission", fontsize=18)
    #ax_dn.set_title("Crosssection of data (at t = 100 ch)", fontsize=24)
    #ax_dn.set_xlabel("x / ch", fontsize=18)
    #ax_dn.set_ylabel("y / ch", fontsize=18)
    #cbar_d = fig.colorbar(im_dn, ax=ax_dn, fraction=0.046, pad=0.04)
    #cbar_d.set_label("Neutron count", fontsize=18)
        ax_1 = fig.add_subplot(grid[1, didx])
        p_d = ax_1.plot(data[:, 100, 36], marker='o', lw=0)
    #ax_d.set_title("2D Map of d$_{110}$", fontsize=24)
    #ax_d1.set_xlabel("x / ch", fontsize=18)
    #ax_d1.set_ylabel("y / ch", fontsize=18)
        ax_1.set_ylim([0, 1.5])
        ax_1.set_ylabel('Transmission', fontsize=18)
        ax_1.set_xlabel('t / ch', fontsize=18)
    #cbar_d = fig.colorbar(im_d, ax=ax_d, fraction=0.046, pad=0.04)
    plt.tight_layout()
    return fig



if __name__ == "__main__":
    transmission_noisy, transmission, transmission_noisy_strd,\
            transmission_strd = get_data()
    figure = plot_expt_transmission(transmission_noisy, transmission,
                                    transmission_noisy_strd, transmission_strd,
                                    cmap="viridis")
    figure.savefig("two_expt_transmission.png", dpi=160, bbox_inches="tight")
    plt.show()

