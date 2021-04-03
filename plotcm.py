#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

params = {'mathtext.default': 'regular', 'axes.linewidth': 1.5}
plt.rcParams.update(params)



def getdata(infile, stepsize):
    m = []
    Cm1 = []
    Cm2 = []
    with open(infile, 'r') as f:
        for line in f:
            values = line.split()
            m.append(float(values[0]))
            Cm1.append(float(values[1]))
            Cm2.append(float(values[2]))
        m = np.array(m)
        Cm1 = np.array(Cm1)/np.prod(stepsize)**2
        Cm2 = np.array(Cm2)/np.prod(stepsize)**2
    return m, Cm1, Cm2


def run():
    stepsize = np.array([0.025, 0.025, 0.025, 0.5])
    infile = "/home/kazu/desktop/200204/fine/hourbyhour/" +\
             "ortho_opt_without_mask/cm.txt"
    m, Cm1, Cm2 = getdata(infile, stepsize)
    alpha = 123.0/5.189990e+05
    plt.plot(m[m > 10.0**5.5]*alpha, Cm1[m > 10.0**5.5], c='gray', ls='dashed',
             label=r'$\alpha$=0, $\Delta_{\omega}$=1.0 meV')
    plt.plot(m[m > 10.0**5.5]*alpha, Cm2[m > 10.0**5.5], c='gray',
             label=r'$\alpha$=0, $\Delta_{\omega}$=1.5 meV')
    infile = "/home/kazu/desktop/200204/fine/hourbyhour/" + \
             "ortho_opt_without_mask/condparam09/cm.txt"
    m, Cm1, Cm2 = getdata(infile, stepsize)
    plt.plot(m[m > 10.0**5.5]*alpha, Cm1[m > 10.0**5.5], c='k', ls='dashed',
             label=r'$\alpha$=0.9, $\Delta_{\omega}$=1.0 meV')
    plt.plot(m[m > 10.0**5.5]*alpha, Cm2[m > 10.0**5.5], c='k',
             label=r'$\alpha$=0.9, $\Delta_{\omega}$=1.5 meV')
    plt.tick_params(top=True, right=True, direction='in', which='both')
    plt.xscale('log')
    plt.xlabel('count in common space')
    plt.ylabel('extrapolated cost function ($meV^{-2} rlu^{-6}$)')
    plt.legend()
    plt.savefig("fig5_costfunc_competetions.pdf")
    plt.show()


run()
