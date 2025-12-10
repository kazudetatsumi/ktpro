#!/usr/bin/evn python
# This script is wrote for simulated data sets of single edge data sets of BI.
# Using 2D patterns of png files, this script generate parameters of single edges energy profiles.
# Kazuyoshi TATSUMI 2025/09/26
import numpy as np
from PIL import Image, ImageFilter
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt


def get_dist(_idx, innerarea, outerarea, area):
    diff = np.linalg.norm(innerarea - _idx, axis=-1)
    edge_inner_pos = innerarea[np.argmin(diff)]
    diff = np.linalg.norm(outerarea - _idx, axis=-1)
    edge_outer_pos = outerarea[np.argmin(diff)]
    dist_inner = np.linalg.norm(_idx - edge_inner_pos)
    dist_outer = np.linalg.norm(_idx - edge_outer_pos)
    return dist_inner, dist_outer


def modify_params(ptn, orgv, params, paramvs,  inner_idx, outer_idx, area_idx, outer_not_area=False):
    innerarea = np.argwhere(ptn == orgv[inner_idx])   # 2
    outerarea = np.argwhere(ptn == orgv[outer_idx])   # 0
    area = np.argwhere(ptn == orgv[area_idx])        # 1
    for _params, _param_values in zip(params, paramvs):
        for _idx in area:
            dist_inner, dist_outer = get_dist(_idx, innerarea, outerarea, area)
            if outer_not_area:
                _params[_idx[0], _idx[1]] = _param_values[outer_idx]*dist_inner/(dist_inner+dist_outer) + _param_values[inner_idx]*dist_outer/(dist_inner+dist_outer)
            else:
                _params[_idx[0], _idx[1]] = _param_values[area_idx]*dist_inner/(dist_inner+dist_outer) + _param_values[inner_idx]*dist_outer/(dist_inner+dist_outer)


def run():
    #ptn1 = np.array(Image.open('denoise_pattern_rev2.png'))[:283, :768]
    #ptn2 = np.array(Image.open('denoise_pattern2_rev.png'))[:283]
    ptn1 = np.array(Image.open('denoise_pattern_rev2.png'))
    ptn2 = np.array(Image.open('denoise_pattern2_rev4.png'))

    #plt.imshow(ptn1, cmap='gray')
    #plt.savefig('org_2d.png')
    #plt.show()
    #plt.plot(ptn1[200])
    #plt.savefig('org_1d.png')
    #plt.show()
    #print(ptn1.shape, ptn2.shape)


    param_values = np.zeros((6, 6))
    params = np.zeros(((6,)+ptn2.shape))
    #param_values[0] = np.array([0.05, 0.2, 0.15, 0.2, 0.])
    #param_values[1] = np.array([0.05, 0.2, 0.15, 0.2, 0.])*(-0.33)
    #param_values[2] = np.array([1.0, 0.5, -0.05, 0.7, 0.])
    #param_values[3] = np.array([1.0, 0.5, -0.05, 0.7, 0.])*(-1.)
    #param_values[0] = np.array([0.05, 0.2, 0.15, 0.2, 0.])*4.
    #param_values[0] = np.array([0.05, 0.1, 0.0001, 0.2, 0.2, 0.])*4.   # A0
    #param_values[1] = param_values[0]*(-0.25)                          # B0
    #param_values[0] = (np.array([0.05, 0.1, 0.0001, 0.2, 0.2, 0.2001])*(-1.) + 0.2001 )*0.1  # A0
    param_values[0] = np.array([0.0150, 0.0100, 0.02000, -0.0200, -0.0200, 0.])*2  # A0
    param_values[1] = np.array([0.05, 0.1, 0.000, 0.2, 0.2, 0.])*4*0.1 # B0
    #param_values[2] = param_values[0]*(-3.0) + 2.                      # Ahkl
    #param_values[2] = param_values[0]*1.5                               # Ahkl
    #param_values[2] = param_values[0]*2.0                               # Ahkl
    #param_values[2] = np.array([0.05, 0.1, 0.0001, 0.2, 0.2, 0.])*8.0    # Ahkl
    param_values[2] = np.array([0.05, 0.1, 0.0001, 0.2, 0.2, 0.])*10.0    # Ahkl
    param_values[3] = param_values[2]*(-0.333)                         # Bhkl
    param_values[4] = np.array([2.036, 2.034, 2.031, 2.028, 2.024, 2.020])
    param_values[5] = np.array([30.0, 25.0, 15.0, 18., 7., 2.])

    edgepos = np.array([[188, 103], [101, 675]])
    det_grad = (edgepos[1, 0] - edgepos[1, 1])/(edgepos[0, 0] - edgepos[1, 0])
    det_ysec = -edgepos[0, 0]*det_grad + edgepos[0, 1]

    orgvalues = np.unique(ptn1)
    print(orgvalues)
    orgvalues = np.unique(ptn2)
    print(orgvalues)
    # [  0  75 100 120 145 255] ptn1
    # [  0  30  55 215 235 255] ptn2
    for sidx in range(6):
        if sidx >= 4:
            ptn = ptn2
        else:
            ptn = ptn1
        for ovidx, orgv in enumerate(np.unique(ptn)):
            params[sidx][ptn == orgv] = param_values[sidx, ovidx]
    plt.imshow(params[0], origin='lower')
    plt.show()
    #cpos = np.array(ptn2.shape)/2-1

    modify_params(ptn2, np.unique(ptn2), params[4:6], param_values[4:6], 3, 5, 0)
    modify_params(ptn2, np.unique(ptn2), params[4:6], param_values[4:6], 1, 3, 2, outer_not_area=True)
    modify_params(ptn2, np.unique(ptn2), params[4:6], param_values[4:6], 4, 2, 1)
    modify_params(ptn1, np.unique(ptn1), params[0:4], param_values[0:4], 1, 5, 2)
    modify_params(ptn1, np.unique(ptn1), params[0:4], param_values[0:4], 0, 2, 1)
    modify_params(ptn1, np.unique(ptn1), params[0:4], param_values[0:4], 3, 0, 4, outer_not_area=True)
    #fig, ax = plt.subplots(6, 2)
    #for pidx, prm in enumerate(params):
    #    ax[pidx, 0].imshow(prm)
    #    ax[pidx, 1].plot(prm[200])
    #plt.show()

    for _params, _param_values in zip(params, param_values):
        for _x in range(_params.shape[0]):
            for _y in range(_params.shape[1]):
                if det_grad*_x + det_ysec < _y:
                    _params[_x, _y] = _param_values[-1]
    y = params[:, 48:-100+8, 100:400]
    plt.imshow(y[0], origin='lower')
    plt.show()
    print(y.shape)
    #print(y.shape)
    y = y[np.newaxis]
    import torch
    y = torch.from_numpy(y)
    upsample = torch.nn.Upsample(size=(72, 192),
                                 mode='bilinear')
    y = upsample(y).reshape((6, 72, 192)).numpy()

    #im = Image.fromarray(params[0])
    #params[0] = np.array(im.filter(ImageFilter.MedianFilter(size = 3)))

    #params[0] = gaussian_filter(params[0], sigma=0.5)
    #plt.imshow(y[0], cmap='gray', origin='lower')
    #plt.savefig('modified_2d.png')
    #plt.show()
    #plt.plot(y[0, 50])
    #plt.savefig('modified_1d.png')
    #plt.show()
    import pickle
    with open('params_scratch_rev3.pkl', 'wb') as f:
        #plt.imshow(y[1])
        #plt.show()
        y[:, y[5] < 3] = 0.
        ## rev2
        #tmp = y[0].copy()  # 1番目のスライスだけコピー
        #plt.imshow(y[1])
        #plt.show()
        #y[0] = y[1]
        #y[1] = tmp*0.1
        #y[0] += 0.3
        #y[:, y[5] < 3] = 0.
        #plt.imshow(y[0])
        #y[0] *= 0.1
        #y[:, y[5] < 3] = 0.
        #plt.imshow(y[0])
        #plt.show()
        pickle.dump(y, f, 4)


run()

#for _idx in np.argwhere(ptn2 == 145):
#    #dist1 = ((_idx[0] - cpos[0])**2 + (_idx[1] - cpos[1])**2)**0.5
#    dist1 = np.linalg.norm(_idx - cpos)
#    for count, __idx in enumerate(np.argwhere(ptn2 == 120)):
#        #dist2 = ((cpos[0] - __idx[0])**2 + (cpos[1] - __idx[1])**2)**0.5 + ((_idx[0] - __idx[0])**2 + (_idx[1] - __idx[1])**2)**0.5
#        dist2 = np.linalg.norm(cpos - __idx)
#        diff = abs(dist1 - dist2)
#        if count == 0:
#            mindiff = diff
#            minidx = __idx
#        elif diff < mindiff:
#            mindiff = diff
#            minidx = __idx
#    edge1 = minidx
#    for count, __idx in enumerate(np.argwhere(ptn2 == 0)):
#        dist2 = np.linalg.norm(cpos - __idx)
#        #dist2 = ((cpos[0] - __idx[0])**2 + (cpos[1] - __idx[1])**2)**0.5 + ((_idx[0] - __idx[0])**2 + (_idx[1] - __idx[1])**2)**0.5
#        diff = abs(dist1 - dist2)
#        if count == 0:
#            mindiff = diff
#            minidx = __idx
#        elif diff < mindiff:
#            mindiff = diff
#            minidx = __idx
#    edge2 = minidx
    #len1 = (_idx[0] - edge1[0])**2 + 




