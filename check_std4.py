#!?usr/bin/env python
import pickle
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
import numpy as np
import os
import socket
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 10
plt.rcParams['mathtext.default'] = 'regular'
mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 9/2*1.3, "axes.labelsize": 11/2*1.3, "xtick.labelsize": 9/2*1.3, "ytick.labelsize": 9/2*1.3,
    "legend.fontsize": 9/2*1.3,
})


def get_tdata():
    tdatafile = 'tdata.pkl'
    if not os.path.exists(tdatafile):
        for i in range(0, 5):
            with open('./seed' + str(i) + '/valtesttot.pkl', 'rb') as f:
                data = pickle.load(f)
                var = pickle.load(f)
            if i == 0:
                tdata = data.squeeze()[np.newaxis]
                tvar = var.squeeze()[np.newaxis]
            else:
                tdata = np.vstack((tdata, data.squeeze()[np.newaxis]))
                tvar = np.vstack((tvar, var.squeeze()[np.newaxis]))
        with open('tdata.pkl', 'wb') as f:
            pickle.dump(tdata.astype('float32'), f, 4)
            pickle.dump(tvar.astype('float32'), f, 4)
    else:
        print('reading ' + tdatafile)
        with open(tdatafile, 'rb') as f:
            tdata = pickle.load(f)
            tvar = pickle.load(f)
    return tdata, tvar


def get_mask(transmission):
    import cv2
    return cv2.GaussianBlur(transmission.sum(axis=-1), (5, 5), 0) > 236


def check_std():
    hostname = socket.gethostname()
    if hostname == "mlfdev61":
        fdir = "/home/kazu/restormer_rev2_lim/bi3d/restormer_conv3d/" +\
             "for_single/train/full/211/true_edge/nll/gau2ch/ktrand/rev4/"
    elif hostname == "mlfdev51":
        fdir = "/data1/kazu/restormer_rev2_lim/bi3d/restormer_conv3d/" +\
             "for_single/train/full/211/true_edge/nll/gau2ch/ktrand/"
    else:
        raise RuntimeError("You should be in mlfdev51 or mlfdev61")
    if not os.path.exists(fdir + 'tmp_data.pkl'):
        tdata, tvar = get_tdata(tdatafile=fdir + 'tdata.pkl')
        M = tdata.shape[0]
        std = (1/M * (tvar.sum(axis=0) + (tdata**2).sum(axis=0)) - np.mean(tdata, axis=0)**2) ** 0.5
        maxstd = np.max(std[:-2], axis=0)
        im = tdata[0, -2].sum(axis=-1)
        transmission = tdata[0, -2] / tdata[0, -1]
        std_sample = std[-2]
        with open(fdir + 'tmp_data.pkl', 'wb') as f:
            pickle.dump(transmission, f, 4)
            pickle.dump(im, f, 4)
            pickle.dump(maxstd, f, 4)
            pickle.dump(std_sample, f, 4)
    else:
        with open(fdir + 'tmp_data.pkl', 'rb') as f:
            transmission = pickle.load(f)
            im = pickle.load(f)
            maxstd = pickle.load(f)
            std_sample = pickle.load(f)

    _mask = get_mask(transmission)
    mask = get_mask(transmission)
    for iy in range(_mask.shape[1]):
        for ix in range(_mask.shape[0]):
            if (mask[ix-3, iy]) & (not mask[ix, iy]):
                _mask[ix, iy] = True

    poslists = [[40, 130], [30, 50], [30, 30], [40, 15], [66, 79]]
    cond = (maxstd < std_sample).sum(axis=-1) + 1
    condin = (maxstd < std_sample).sum(axis=-1)
    cond[mask] = 2
    condin[_mask] = 0.
    maxpos = np.unravel_index(np.argmax(cond, axis=None), cond.shape)
    maxposin = np.unravel_index(np.argmax(condin, axis=None), cond.shape)

    # ===== layout in inch ====
    axes_width = 3.0*0.6     # width of map or lineplots
    map_height = 2.2*0.6     # height of map in the 1st row
    lineplot_height = 1.2*3*0.6   # total hiehgt of three lineplots in the 2nd row
    left_margin = 0.6*0.6    # left vacant space
    right_margin = 0.6*0.6   # right vacant space in the righter side of the cbar
    top_margin = 0.3*0.6     # top vacant space
    bottom_margin = 0.4*0.6  # bottom vacant space
    # ===== colorbar =====
    cbar_width = 0.05     # width of cbar
    cbar_pad = 0.02       # horizontal distacne bewteen map and colorbar
    # ===== vacant space between the rows =====
    hspace = 0.0
    # ===== figure ==========================
    fig_height = top_margin + map_height + lineplot_height + bottom_margin
    fig_width = left_margin + axes_width + cbar_pad + cbar_width + right_margin
    print(fig_width, fig_height)
    fig = plt.figure(figsize=(fig_width, fig_height))
    # ====== grid spec (2x1) =========================
    #   1st row: 2Dmap, 2nd row: lineplots which will be divided to three axes
    #   "cbar_pad + cbar_width" is kept in the right end  for colorbar
    # ===============================
    height_ratios = [map_height, lineplot_height]
    gs = GridSpec(
        2, 1, figure=fig, height_ratios=height_ratios,
        left=left_margin/fig_width,
        right=1 - (right_margin + cbar_pad + cbar_width)/fig_width,
        bottom=bottom_margin/fig_height,
        top=1 - top_margin/fig_height,
        hspace=hspace
    )

    # ===== 1st row of 2D map=====
    ax_map = fig.add_subplot(gs[0, 0])
    cim = ax_map.imshow(cond, norm=colors.LogNorm(vmin=cond.min(),
                        vmax=cond.max()), cmap='Grays', origin='lower')

    # Annotation
    for pidx, pos in enumerate([maxpos, maxposin, poslists[0]]):
        ax_map.annotate('#'+str(pidx+1), xy=(pos[1], pos[0]),
                        xytext=(pos[1]-20, pos[0]-20), textcoords='data',
                        color='k',
                        arrowprops=dict(arrowstyle="->", color='k'))

    ax_map.set_ylabel('y / ch')
    ax_map.set_xlabel('x / ch')
    ax_map.tick_params(direction='in')

    # ===== put cbar in the right side of the map strictly by add_axes =====
    map_pos = ax_map.get_position()  # figure coordination
    cbar_width_fig = cbar_width / fig_width
    cbar_pad_fig = cbar_pad / fig_width
    cax = fig.add_axes([
        map_pos.x1 + cbar_pad_fig,  # left
        map_pos.y0,                 # bottom = y0 of map
        cbar_width_fig,             # width
        map_pos.height              # height = height of map
    ])
    cbar = fig.colorbar(cim, cax=cax)
    cbar.set_label("# of violated voxels + 1")
    cax.tick_params(direction='out', which='both')

    # ===== 2nd row: three lineplots without no vacant space between them =====
    #area_pos = fig.add_subplot(gs[1, 0]).get_position()
    area_pos = gs[1, 0].get_position(fig)

    # equally divided into three and no vacant space between them.
    h_each = area_pos.height / 3.0
    left = area_pos.x0
    width = area_pos.width
    bottom4 = area_pos.y0
    bottom3 = bottom4 + h_each
    bottom2 = bottom3 + h_each
    ax4 = fig.add_axes([left, bottom4, width, h_each])
    ax3 = fig.add_axes([left, bottom3, width, h_each], sharex=ax4)
    ax2 = fig.add_axes([left, bottom2, width, h_each], sharex=ax4)
    tof = np.arange(maxstd.shape[-1])*20 + 23000
    for axid, (axi, pos) in enumerate(zip([ax2, ax3, ax4],
                                          [maxpos, maxposin, poslists[0]])):
        axi.plot(tof, maxstd[pos[0], pos[1]], label='Ref.', c='k', ls='--')
        axi.plot(tof, std_sample[pos[0], pos[1]], label='Expt.', c='k')
        axi.tick_params(direction='in', top=True, right=True)
        axi.set_ylim([0., 6])
        axi.set_xlim([tof[0], tof[-1]])
        axi.legend(loc='upper right', frameon=True)
        axi.set_ylabel('std')
        axi.text(2, 5, "#"+str(axid+1))
    for axi in [ax2, ax3]:
        plt.setp(axi.get_xticklabels(), visible=False)
        axi.tick_params(labelbottom=False)
    ax4.set_xlabel(r'TOF / $\mu$s')
    for axi in [ax3, ax4]: # 下側の軸の「一番上の目盛り」を消す場合
        axi.yaxis.get_major_ticks()[-1].label1.set_visible(False)
    plt.savefig('fig_bi_std.eps')
    plt.show()


if __name__ == '__main__':
    check_std()
