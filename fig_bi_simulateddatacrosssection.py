import numpy as np
import copy
import matplotlib.pyplot as plt
import matplotlib as mpl

#import matplotlib as mpl
plt.rcParams["font.family"] = "Arial"   # 使用するフォント
plt.rcParams["font.size"] = 9
plt.rcParams['mathtext.default'] = 'regular'

# フォント・EPS設定
mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 9, "axes.labelsize": 11, "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'


def get_data():
    import pickle
    with open('/home/kazu/desktop/240424/connect2d/params_scratch_rev4.pkl',
              'rb') as f:
        y = pickle.load(f)
    with open('/home/kazu/desktop/240424/connect2d/test_single_edge/" + \
              "bi3d_scratch_rev4_211_phantom_local.pkl',
              'rb') as f:
        sample = pickle.load(f)
    return y, sample


def plot_three_maps_with_colorbars_and_1d_profiles(
    y,
    cmap: str = "viridis",
) -> plt.Figure:
    fig = plt.figure(figsize=(8, 12))
    grid = fig.add_gridspec(
        nrows=y.shape[0],
        ncols=2,
        height_ratios=np.ones((y.shape[0])),
        width_ratios=[0.63, 0.37],
        hspace=0.5,
        wspace=0.12*2,
    )
    variables = ['Neutron Count', r'$\mathrm{a_0}$', r'$\mathrm{b_0}$',
                 r'$\mathrm{a_{110}}$', r'$\mathrm{b_{110}}$',
                 r'$\mathrm{d_{110}}$', r'$\mathrm{{\sigma}_{0}}$']
    units = ['',  '', r' / $\mu\mathrm{s^{-1}}$', '',
             r' / $\mu\mathrm{s^{-1}}$', r' / $\AA$',
             r' / $\mu\mathrm{s^{-1}}$', ]
    for yidx, (_y, ystr, yunit) in enumerate(zip(y, variables, units)):

        __y = copy.deepcopy(_y)
        __y[__y == 0.] = np.nan
        cmap = copy.copy(plt.get_cmap('gray_r'))
        cmap.set_bad(color='0.9')
        ax = fig.add_subplot(grid[yidx, 0])
        if ystr == r'$\mathrm{d_{110}}$':
            im = ax.imshow(
                __y.T,
                origin="lower",
                cmap=cmap,
                aspect=1,
                interpolation="nearest",
                vmin=np.unique(np.sort(__y))[1]
            )
        else:
            im = ax.imshow(
                __y.T,
                origin="lower",
                cmap=cmap,
                aspect=1,
                interpolation="nearest",
            )
        ax.set_xlabel("x / ch", labelpad=0)
        ax.set_ylabel("y / ch", labelpad=0)
        if yidx == 0:
            ax.text(2, 2, f"{ystr}", c='w')
        else:
            ax.text(2, 2, f"{ystr}", c='k')
        ax.axhline(40, color='k', linestyle='-', alpha=0.8, lw=2.0)
        ax.axhline(40, color='w', linestyle='-', lw=1.0)
        ax.tick_params(direction='in', length=3, pad=1.4)
        ax.set_xticks([0, 50, 100, 150])

        cbar_d = fig.colorbar(im, ax=ax, pad=0.008)
        cbar_d.set_label(f"{ystr}" + f"{yunit}", labelpad=1)
        cbar_d.ax.tick_params(pad=2, direction='out')

        ax1 = fig.add_subplot(grid[yidx, 1])
        ax1.plot(_y.T[40], c='k')
        ax1.set_xlim([0, y.shape[1]])
        ax1.tick_params(direction="in")
        ax1.set_xticks([0, 50, 100, 150])
        ax1.set_ylabel(f"{ystr}"+f"{yunit}", labelpad=0)
        ax1.set_xlabel("x / ch", labelpad=0)
        ax1.tick_params(length=3, pad=2)
        if ystr == r'$\mathrm{d_{110}}$':
            ax1.set_ylim([np.unique(np.sort(_y))[1],
                          np.unique(np.sort(_y))[-1]])
    return fig


if __name__ == "__main__":
    y, sample = get_data()
    y = y.transpose((0, 2, 1))
    D = sample[:, 2:-2, 2:-2].transpose((1, 2, 0))[:, :, 100]
    Dwy = np.vstack((D[np.newaxis], y))

    figure = plot_three_maps_with_colorbars_and_1d_profiles(Dwy)
    plt.show()
