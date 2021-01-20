#!/usr/bin/env python
import numpy as np
P = 0
S = 2
FH = 2
FW = 2
N = 1
C = 1
H = 12
W = 14
OH = (H + 2*P - FH)//S +1
OW = (W + 2*P - FW)//S +1
col = np.ones((N, C, FH, FW, OH, OW))
img = np.zeros((N, C, H+2*P + S -1, W+2*P + S -1))


for y in range(FH):
    y_max = y + S*OH
    for x in range(FW):
        x_max = x + S*OW
        img[:, :, y:y_max:S, x:x_max:S] += col[:, :, y, x, :, :]


print(img)
print(img[:,:,P:P+H, P:P+W])
print(H+2*P-FH)
