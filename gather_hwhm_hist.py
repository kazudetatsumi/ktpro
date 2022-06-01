#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as so
#import os

params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.linewidth'] = 1.2 


def run():
    prefix = "/home/kazu/desktop/210108/Tatsumi"
    runs = ["6202", "6205", "6203", "6206", "6207"]
    temps = np.array([303, 288, 275, 263, 253])
    logfiles = [prefix+"/winparam_exam_6202/160_1_0000025io_Boxcar/n8000/" +
                "std-0000025io-160-1-Boxcar-2lore.log"]
    for run, temp in zip(runs[1:], temps[1:]):
        logfile = prefix + "/winparam_exam_" + run +\
                  "/160_1_0000001io_Boxcar/n200000/" + \
                  "std-0000001io-160-1-Boxcar-2lore.log"
        logfiles.append(logfile)
        # print(os.path.isfile(logfile))
    hwhm = get_hwhms(logfiles)
    print(hwhm)
    for hidx in range(0, 2):
        [k0, k1], R2 = optimize(1./temps, np.log(hwhm[:, hidx]),
                                variables=[0.0, -1.0])
        x = np.linspace(np.min(1./temps), np.max(1./temps), 100)
        plt.scatter(1./temps, np.log(hwhm[:, hidx]),
                    label=r'$\Gamma$'+str(hidx))
        plt.plot(x, k0 + k1*x)
        plt.text(0.00355, -5.2-hidx*3.3,  "$R^2$={:.3f}".format(R2))
        plt.tick_params(direction='in', right=True, top=True, labelbottom=True,
                        width=1.2)
        plt.xlabel('1/T [$K^{-1}$]')
        plt.ylabel('log($\Gamma$)')

        plt.legend() 
    plt.show()


def res(coeffs, x, t):
    [k0, k1] = coeffs
    y = k0 + k1*x
    return t-y


def optimize(x, t, variables=[0.0, -1.0]):
    out = so.least_squares(res, variables, args=(x, t))
    R2 = 1 - 2. * out.cost  /(np.var(t)*t.shape[0])
    print(R2)
    return out.x, R2


def get_hwhms(logfiles):
    hwhm = np.zeros((len(logfiles), 2))
    for ilog, logfile in enumerate(logfiles):
        with open(logfile, 'r') as f:
            lines = f.readlines()
            for il, line in enumerate(lines):
                if 'cov' in line:
                    hwhm[ilog, 0] = lines[il-1].split( )[1]
                    hwhm[ilog, 1] = lines[il-1].split( )[3]
                    for hidx in range(0, 2):
                        if hwhm[ilog, hidx] < 0.:
                            hwhm[ilog, hidx] = hwhm[ilog, hidx]*(-1.)
                    if hwhm[ilog, 0] < hwhm[ilog, 1]:
                        tmpvalue = hwhm[ilog, 0]
                        hwhm[ilog, 0] = hwhm[ilog, 1]
                        hwhm[ilog, 1] = tmpvalue
    return(hwhm)




run()
    

