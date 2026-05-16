#!/user/bin/env python
import matplotlib.pyplot as plt
import numpy as np

params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.linewidth'] = 1.2


def plotdata(spectrafile):
    data = np.loadtxt(open(spectrafile, 'rb'), skiprows=1).T

    plt.figure(figsize=(8, 6))

    plt.scatter(data[0], data[2], label='target', marker='x', c='k')
    plt.plot(data[0], data[3], label='ML', c='k')
    plt.legend()
    plt.xlabel(u'TOF ${\mu}s$')
    plt.ylabel(u'INTENSITY')
    plt.tick_params(direction='in', right=True, top=True, labelbottom=True,
                    width=1.2)
    plt.show()




def run():
    spectrafile = "/home/kazu/desktop/220507/raim/gpr/output/ca_pos_ta/" +\
                  "final_ta.txt"
    plotdata(spectrafile)


run()
