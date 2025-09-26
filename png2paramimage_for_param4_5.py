#!/usr/bin/evn python
# This script is wrote for simulated data sets of single edge data sets of BI.
# Using 2D patterns of png files, this script generate parameters of single edges energy profiles.
# Kazuyoshi TATSUMI 2025/09/26
import numpy as np
from PIL import Image, ImageFilter
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt

ptn1 = np.array(Image.open('denoise_pattern_rev.png'))
ptn2 = np.array(Image.open('denoise_pattern2_rev.png'))

plt.imshow(ptn2, cmap='gray')
plt.savefig('org_2d.png')
plt.show()
plt.plot(ptn2[200])
plt.savefig('org_1d.png')
plt.show()

centerorgv = 120
a0_center = 0.1
a0_nearcenter = 0.2
a0_valley = 0.05
a0_edge = 0.2

param_values = np.zeros((6, 5))
params = np.zeros(((6,)+ptn2.shape))
param_values[0] = np.array([0.05, 0.2, 0.15, 0.2, 0.])
param_values[1] = np.array([0.05, 0.2, 0.15, 0.2, 0.])*(-0.33)
param_values[2] = np.array([1.0, 0.5, -0.05, 0.7, 0.])
param_values[3] = np.array([1.0, 0.5, -0.05, 0.7, 0.])*(-1.)
param_values[4] = np.array([2.028, 2.036, 2.032, 2.024, 2.020])
param_values[5] = np.array([30.0, 25.0, 15.0, 7., 2.])

edgepos = np.array([[188, 103], [101, 675]])
det_grad = (edgepos[1, 0] - edgepos[1, 1])/(edgepos[0, 0] - edgepos[1, 0])
det_ysec = -edgepos[0, 0]*det_grad + edgepos[0, 1]

orgvalues = np.unique(ptn2)
print(orgvalues)
# [  0 100 120 145 255]
# [  0  75 100 235 255]
for sidx in range(6):
    for ovidx, orgv in enumerate(orgvalues):
        params[sidx][ptn2 == orgv] = param_values[sidx, ovidx]
cpos = np.array(ptn2.shape)/2-1

innerarea = np.argwhere(ptn2 == orgvalues[1])   # 2
outerarea = np.argwhere(ptn2 == orgvalues[4])   # 0
area = np.argwhere(ptn2 == orgvalues[0])        # 1
dist = np.linalg.norm(area - cpos, axis=-1)
for _params, _param_values in zip(params[5:6], param_values[5:6]):
    for _idx, _dist in zip(area, dist):
        diff = np.linalg.norm(innerarea - _idx, axis=-1)
        edge_inner_pos = innerarea[np.argmin(diff)]
        diff = np.linalg.norm(outerarea - _idx, axis=-1)
        edge_outer_pos = outerarea[np.argmin(diff)]
        edge_dist = np.linalg.norm(edge_outer_pos - edge_inner_pos)
        dist_inner = np.linalg.norm(_idx - edge_inner_pos)
        dist_outer = np.linalg.norm(_idx - edge_outer_pos)
        _params[_idx[0], _idx[1]] = _param_values[0]*dist_inner/(dist_inner+dist_outer) + _param_values[1]*dist_outer/(dist_inner+dist_outer)
innerarea = np.argwhere(ptn2 == orgvalues[2])
outerarea = np.argwhere(ptn2 == orgvalues[0])
area = np.argwhere(ptn2 == orgvalues[1])
dist = np.linalg.norm(area - cpos, axis=-1)
for _params, _param_values in zip(params[5:6], param_values[5:6]):
    for _idx, _dist in zip(area, dist):
        diff = np.linalg.norm(innerarea - _idx, axis=-1)
        edge_inner_pos = innerarea[np.argmin(diff)]
        diff = np.linalg.norm(outerarea - _idx, axis=-1)
        edge_outer_pos = outerarea[np.argmin(diff)]
        edge_dist = np.linalg.norm(edge_outer_pos - edge_inner_pos)
        dist_inner = np.linalg.norm(_idx - edge_inner_pos)
        dist_outer = np.linalg.norm(_idx - edge_outer_pos)
        _params[_idx[0], _idx[1]] = _param_values[1] * dist_inner / (dist_inner + dist_outer) + _param_values[2] * dist_outer / (dist_inner + dist_outer)
innerarea = np.argwhere(ptn2 == orgvalues[3])
outerarea = np.argwhere(ptn2 == orgvalues[1])
area = np.argwhere(ptn2 == orgvalues[2])
dist = np.linalg.norm(area - cpos, axis=-1)
for _params, _param_values in zip(params[5:6], param_values[5:6]):
    for _idx, _dist in zip(area, dist):
        diff = np.linalg.norm(innerarea - _idx, axis=-1)
        edge_inner_pos = innerarea[np.argmin(diff)]
        diff = np.linalg.norm(outerarea - _idx, axis=-1)
        edge_outer_pos = outerarea[np.argmin(diff)]
        edge_dist = np.linalg.norm(edge_outer_pos - edge_inner_pos)
        dist_inner = np.linalg.norm(_idx - edge_inner_pos)
        dist_outer = np.linalg.norm(_idx - edge_outer_pos)
        _params[_idx[0], _idx[1]] = _param_values[2] * dist_inner / (dist_inner + dist_outer) + _param_values[3] * dist_outer / (dist_inner + dist_outer)
for _params, _param_values in zip(params[5:6], param_values[5:6]):
    for _x in range(_params.shape[0]):
        for _y in range(_params.shape[1]):
            if det_grad*_x + det_ysec < _y:
                _params[_x, _y] = _param_values[-1]
y = params[5, 40:-100, 100:400]
print(y.shape)
#print(y.shape)
y = y[np.newaxis, np.newaxis]
import torch
y = torch.from_numpy(y)
upsample = torch.nn.Upsample(size=(72, 192),
                             mode='bilinear')
y = upsample(y).reshape((72, 192)).numpy()

im = Image.fromarray(params[0])
#params[0] = np.array(im.filter(ImageFilter.MedianFilter(size = 3)))

#params[0] = gaussian_filter(params[0], sigma=0.5)
plt.imshow(y, cmap='gray', origin='lower')
plt.savefig('modified_2d.png')
plt.show()
plt.plot(y[50])
plt.savefig('modified_1d.png')
plt.show()
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




