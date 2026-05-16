import numpy as np
import matplotlib.pyplot as plt

#import matplotlib as mpl
plt.rcParams["font.family"] = "Arial"   # 使用するフォント
plt.rcParams["font.size"] = 9
plt.rcParams['mathtext.default'] = 'regular'

## Make Arial the default font for all text
#mpl.rcParams['font.family'] = 'Arial'         # primary family
#mpl.rcParams['font.sans-serif'] = ['Arial']   # optional: ensure sans-serif points to Arial
#
## Optional: math text style (see section 4)
#mpl.rcParams['mathtext.fontset'] = 'dejavusans'  # keeps math readable with a sans style
#

def get_data():
    import pickle
    with open('params_scratch_rev4.pkl', 'rb') as f:
        y = pickle.load(f)
    with open('./test_single_edge/bi3d_scratch_rev4_211_phantom_local.pkl', 'rb') as f:
        sample = pickle.load(f)
    return y, sample


def plot_three_maps_with_colorbars(
    D: np.ndarray,
    P1: np.ndarray,
    P2: np.ndarray,
    cmap: str = "viridis",
) -> plt.Figure:
    # Figure and grid (3 rows, 1 column)
    fig = plt.figure(figsize=(8, 13))
    grid = fig.add_gridspec(
        nrows=3,
        ncols=1,
        height_ratios=[1.0, 1.0, 1.0],
        hspace=0.35,
    )

    # Top: D map + colorbar
    ax_d = fig.add_subplot(grid[0, 0])
    im_d = ax_d.imshow(
            #d_sum_z.T,
        D.T,
        origin="lower",
        cmap=cmap,
        aspect="auto",
        interpolation="nearest",
    )
    ax_d.set_title("Crosssection of data (at t = 100 ch)", fontsize=24)
    ax_d.set_xlabel("x / ch", fontsize=18)
    ax_d.set_ylabel("y / ch", fontsize=18)
    cbar_d = fig.colorbar(im_d, ax=ax_d, fraction=0.046, pad=0.04)
    cbar_d.set_label("Neutron count", fontsize=18)

    # Middle: P1 map + colorbar (share axes with D)
    ax_p1 = fig.add_subplot(grid[1, 0], sharex=ax_d, sharey=ax_d)
    im_p1 = ax_p1.imshow(
        P1.T,
        origin="lower",
        cmap=cmap,
        aspect="auto",
        interpolation="nearest",
        vmin=2.02,
        vmax=2.04
    )
    ax_p1.set_title("2D Map of d$_{110}$", fontsize=24)
    ax_p1.set_xlabel("x / ch", fontsize=18)
    ax_p1.set_ylabel("y / ch", fontsize=18)
    cbar_p1 = fig.colorbar(im_p1, ax=ax_p1, fraction=0.046, pad=0.04)
    cbar_p1.set_label("d$_{110}$ / $\\mathrm{\\AA}$", fontsize=18)

    # Bottom: P2 map + colorbar (share axes with D)
    ax_p2 = fig.add_subplot(grid[2, 0], sharex=ax_d, sharey=ax_d)
    im_p2 = ax_p2.imshow(
        P2.T,
        origin="lower",
        cmap=cmap,
        aspect="auto",
        interpolation="nearest",
    )
    ax_p2.set_title("2D Map of a$_0$", fontsize=24)
    ax_p2.set_xlabel("x / ch", fontsize=18)
    ax_p2.set_ylabel("y / ch", fontsize=18)
    cbar_p2 = fig.colorbar(im_p2, ax=ax_p2, fraction=0.046, pad=0.04)
    cbar_p2.set_label( "a$_0$", fontsize=18)

    plt.tight_layout()
    return fig


def plot_three_maps_with_colorbars_and_1d_profiles(
    y,
    cmap: str = "viridis",
) -> plt.Figure:
    # Figure and grid (3 rows, 1 column)
    fig = plt.figure(figsize=(8, 12))
    #fig.suptitle('Distributions of neutron counts and parameters in the simulated data point', fontsize=12)
    grid = fig.add_gridspec(
        nrows=y.shape[0],
        ncols=2,
        height_ratios=np.ones((y.shape[0])),
        width_ratios=[0.63, 0.37],
        hspace=0.5,
        wspace=0.12,
    )
    variables = ['Neutron Count', r'$\mathrm{a_0}$', r'$\mathrm{b_0}$', r'$\mathrm{a_{110}}$',r'$\mathrm{b_{110}}$',r'$\mathrm{d_{110}}$',r'$\mathrm{{\sigma}_{0}}$']
    units = ['',  '', r' / $\mu\mathrm{s^{-1}}$', '', r' / $\mu\mathrm{s^{-1}}$', r' / $\AA$', r' / $\mu\mathrm{s^{-1}}$', ]
    for yidx, (_y, ystr, yunit) in enumerate(zip(y, variables, units)):
            #for yidx, (_y, ystr) in enumerate(zip(y, ['Neutron Count', r'$\mathrm{a_0}$',
            #r'$\mathrm{b_0 \ /\ \mu}\mathrm{s^{-1}}$', r'$\mathrm{a_{110}}$',
            #r'$\mathrm{b_{110}\ /\ \mu}\mathrm{s^{-1}}$',
            #r'$\mathrm{d_{110}\ /\ \AA}$',
            #r'$\mathrm{{\sigma}_{0}\ /\ \mu}\mathrm{s^{-1}}$'])):
        ax = fig.add_subplot(grid[yidx, 0])
        if ystr == r'$\mathrm{d_{110}}$':
            im = ax.imshow(
                _y.T,
                origin="lower",
                cmap=cmap,
                aspect=1,
                interpolation="nearest",
                vmin=np.unique(np.sort(_y))[1]
            )
        else:
            im = ax.imshow(
                _y.T,
                origin="lower",
                cmap=cmap,
                aspect=1,
                interpolation="nearest",
            )
        #ax.set_title(f"{ystr}", fontsize=12)
        ax.set_xlabel("x / ch", labelpad=0)
        ax.set_ylabel("y / ch", labelpad=0)
        ax.text(2, 2, f"{ystr}", c='w')
        ax.axhline(40, color='r', linestyle='--', lw=1)
        ax.tick_params(length=3, labelsize=8, pad=1.4)
        ax.set_xticks([0, 50, 100, 150])


        cbar_d = fig.colorbar(im, ax=ax, pad=0.005)
        cbar_d.set_label(f"{ystr}" + f"{yunit}", fontsize=8, labelpad=1)
        cbar_d.ax.tick_params(labelsize=8, pad=2, direction='in')

        ax1 = fig.add_subplot(grid[yidx, 1])
        ax1.plot(_y.T[40])
        ax1.set_xlim([0, y.shape[1]])
        ax1.tick_params(direction = "in", labelsize=8)
        ax1.set_xticks([0, 50, 100, 150])
        ax1.set_ylabel(f"{ystr}"+f"{yunit}", labelpad=0)
        ax1.set_xlabel("x / ch", labelpad=0)
        ax1.tick_params(length=3, labelsize=8, pad=2)
        #if yidx == y.shape[0]-1:
        #    ax1.set_xlabel("x / ch")
        #    ax1.tick_params(direction = "in")
        #else:
        #    ax1.tick_params(labelbottom=False)
        #    ax1.tick_params(direction = "in")
        #leg = ax1.legend(frameon=False, handlelength=1.2,
        #        title=f'at y=40', alignment='left')
        #leg.set_title(f'at y=40')

        if ystr == r'$\mathrm{d_{110}}$':
            ax1.set_ylim([np.unique(np.sort(_y))[1], np.unique(np.sort(_y))[-1]])
    return fig


if __name__ == "__main__":
    y, sample = get_data()
    y = y.transpose((0, 2, 1))
    D = sample[:, 2:-2, 2:-2].transpose((1, 2, 0))[:,:,100]
    print(y.shape)
    print(D.shape)
    Dwy = np.vstack((D[np.newaxis], y))
    print(Dwy.shape)
    Nx = y.shape[0]
    Ny = y.shape[1]
    Nz = y.shape[2]

    P1 = y[-2]
    P2 = y[0]

    #figure = plot_three_maps_with_colorbars(D, P1, P2, cmap="viridis")
    figure = plot_three_maps_with_colorbars_and_1d_profiles(Dwy)
    #figure.savefig("three_maps_with_colorbars.png", dpi=160, bbox_inches="tight")
    #figure.savefig("six_parameters_maps_with_colorbars_and_1D_profiles.png", dpi=160, bbox_inches="tight")
    figure.savefig("bi_simulated_params_wcb.eps")
    plt.show()

