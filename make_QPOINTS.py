#!/usr/bin/env python
import numpy as np
numqp = 10
qp = np.zeros((numqp, 3))

for i in range(0, numqp):
    qp[i, 0] = 0
    qp[i, 1] = 0
    qp[i, 2] = float(i) / numqp


with open('./qpoiintfile_tmp', mode='w') as w:
    w.write("%d \n" % numqp)
    w.write('\n'.join(["%20.17f %20.17f %20.17f" % tuple(q) for q in qp]))
