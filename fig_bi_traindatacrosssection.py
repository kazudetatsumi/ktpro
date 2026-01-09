#!/usr/bin/env python
import pickle
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Arial"   # 使用するフォント
plt.rcParams["font.size"] = 16
plt.rcParams['mathtext.default'] = 'regular'


def synthesize_bi3ddata():
    with open('bi3d_testbcc_simudata_rev2_lim_single_resize_full_211_true_edgev_ktrand_wosingle.pkl', 'rb') as f:
        sample_saved = pickle.load(f).squeeze().transpose((0, 3, 1, 2))
    _shape = sample_saved.shape
    x = np.arange(_shape[1])*20.+23000
    py = 57
    px = 44
    snum = 4
    fig, ax = plt.subplots(2, figsize=(7, 7))
    _sample = np.random.poisson(sample_saved[snum])
    ax[0].set_title('Simulated sample projected onto (x, y)')
    ax[0].imshow(_sample.sum(axis=0), origin='lower')
    ax[0].set_xlabel('x / ch')
    ax[0].set_ylabel('y / ch')
    ax[0].annotate(f'({py}, {px})', xy=(py, px),
                   xytext=(py+5, px+5), textcoords='data',
                   color='white',
                   arrowprops=dict(arrowstyle="->", color='white'))
    ax[1].set_title(f'Simulated sample intensity at x={py} ch, y={px} ch')
    ax[1].plot(x, _sample[:, px, py])
    ax[1].set_xlabel(r'TOF / $\mu$s')
    ax[1].set_ylabel('Neutron Count')
    plt.tight_layout()
    plt.savefig('bi_simulated_profile.eps')
    plt.show()

synthesize_bi3ddata()
