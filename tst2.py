#!/usr/bin/env python
import matplotlib.pyplot as plt
import math
import matplotlib
import numpy as np



N=128
f = 8
omg = [i*f*2*math.pi/N for i in range(N)]
sig = [math.sin(i) for i in omg]
#cind = [0.1 * x for x in range(N)]
cind = np.random.rand(N)*np.array(sig)

fig,(ax) = plt.subplots(1,1)
colormap = plt.get_cmap('jet')
#norm = matplotlib.colors.Normalize(vmin=min(sig), vmax=max(sig))
#c = colormap(norm(sig))
norm = matplotlib.colors.Normalize(vmin=min(cind), vmax=max(cind))
c = colormap(norm(cind))

for j in range(len(sig)-1):
        ax.plot(omg[j:j+2], sig[j:j+2], color=c[j+1])

plt.show()
