#!/usr/bin/env python
import sys
import os
import numpy as np
sys.path.append(os.pardir)
from common.util import im2col


X = np.random.rand(30, 10, 100, 200)
W = np.random.rand(7, 10, 2, 2)
N, C, H, w = X.shape
FN, C, FH, FW = W.shape
S = 1
P = 0

OH = int(1 + (H  + 2*P - FH) // S)
OW = int(1 + (w  + 2*P - FW) // S)
print(OH, OW)


col = im2col(X, FH, FW, S, P)
col_W = W.reshape(FN, -1).T


out = np.dot(col, col_W)
out = out.reshape(N, OH, OW, FN).transpose(0, 3, 1 ,2)

print(out.shape)



