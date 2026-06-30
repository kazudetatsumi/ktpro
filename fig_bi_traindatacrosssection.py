#!/usr/bin/env python
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import ScalarFormatter

plt.rcParams["font.family"] = "Arial"   # 使用するフォント
plt.rcParams["font.size"] = 16
plt.rcParams['mathtext.default'] = 'regular'

mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 9, "axes.labelsize": 11, "xtick.labelsize": 9, "ytick.labelsize": 9,
    "legend.fontsize": 9,
})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'


def synthesize_bi3ddata():
    fdir = "/home/kazu/desktop/240424/connect2d/sampled/"
    pklfile = "bi3d_testbcc_simudata_rev2_lim_single_resize_full_211_true_edgev_ktrand.pkl"
    with open(fdir+pklfile, 'rb') as f:
        sample_saved = pickle.load(f).squeeze().transpose((0, 3, 1, 2))
    _shape = sample_saved.shape
    x = np.arange(_shape[1])*20.+23000
    py = 57
    px = 44
    py = 50
    px = 50
    snum = 4
    fig, ax = plt.subplots(2, figsize=(4, 4))
    _sample = np.random.poisson(sample_saved[snum])
    ax[0].set_title('Neutron counts of training data point projected onto (x, y)')
    im = ax[0].imshow(_sample.sum(axis=0), origin='lower', cmap='gray')
    #im = ax[0].imshow(_sample.sum(axis=0), origin='lower')
    ax[0].set_xlabel('x / ch')
    ax[0].set_ylabel('y / ch')
    ax[0].tick_params(direction='in', color='w')
    ax[0].annotate('', xy=(py, px),
                   xytext=(py+7, px+7), textcoords='data',
                   color='white',
                   arrowprops=dict(arrowstyle="->", color='white'))
    #cbar = fig.colorbar(im)
    #fmt = ScalarFormatter(useMathText=True)
    #fmt.set_powerlimits((-3, 3))  # 10^n を使う境界。必要に応じて調整
    #cbar.formatter = fmt
    #cbar.ax.yaxis.set_offset_position('right')
    #cbar.update_ticks()

    ax[1].set_title(f'TOF spectrum of training data point at x={py} ch, y={px} ch')
    ax[1].plot(x, _sample[:, px, py], marker='o', ms=2, lw=0, c='k')
    ax[1].plot(x, sample_saved[snum][:, px, py], ls="-", c='gray')

    ax[1].set_xlabel(r'TOF / $\mu$s')
    ax[1].set_ylabel('Neutron Count')
    ax[1].tick_params(direction='in')
    ax[1].set_xlim([x[0], x[-1]])
    plt.tight_layout()
    map_pos = ax[0].get_position()
    cax = fig.add_axes([
        map_pos.x1 + 0.005,  # left
        map_pos.y0,                 # bottom（マップと一致）
        0.010,             # width
        map_pos.height              # height（マップと完全一致）
    ])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label('Neutron Count', labelpad=0.1)
    cax.tick_params(direction='out', labelsize=8, length=2)
    cbar.locator = MaxNLocator(nbins=4)
    cbar.update_ticks()
    plt.savefig('fig_bi_trainingdata_.eps')
    plt.show()


synthesize_bi3ddata()
