#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py
import ctypes
lib = ctypes.CDLL("./histfort_.so")


def get2ddata(f, xi, xf, yi, yf):
    data = np.genfromtxt(f,  delimiter=',', dtype=None)
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    dx = 0.005
    dy = 0.1
    xlin = np.arange(min(x), max(x)+dx, dx)
    nx = xlin.shape[0]
    ylin = np.arange(min(y), max(y)+dy, dy)
    ny = ylin.shape[0]
    karr = np.zeros((nx, ny))
    karr2 = np.zeros((nx, ny))

    for _x, _y, _z in zip(x, y, z):
        xx = np.where(abs(xlin - _x) < 0.0000001)
        yy = np.where(abs(ylin - _y) < 0.0000001)
        karr[xx, yy] = _z + 0.00000001

    condition = karr > 0.0000000001
    karrnonzero = np.extract(condition, karr)

    for _x, _y, _z in zip(x, y, z):
        xx = np.where(abs(xlin - _x) < 0.0000001)
        yy = np.where(abs(ylin - _y) < 0.0000001)
        karr2[xx, yy] = _z

    return karr2[xi:xf, yi:yf], condition[xi:xf, yi:yf]


def get2ddata_for_commandline(f, xi, xf, yi, yf):
    data = np.genfromtxt(f,  delimiter=',', dtype=None)
    x = np.extract(data[:, 2] < 1e+100, data[:, 0])
    y = np.extract(data[:, 2] < 1e+100, data[:, 1])
    z = np.extract(data[:, 2] < 1e+100, data[:, 2])
    dx = 0.005
    dy = 0.1
    xlin = np.arange(min(x), max(x)+dx, dx)
    nx = xlin.shape[0]
    ylin = np.arange(min(y), max(y)+dy, dy)
    ny = ylin.shape[0]
    karr = np.zeros((nx, ny)) 
    karr2 = np.zeros((nx, ny)) 

    for _x, _y, _z in zip(x, y, z):
        xx = np.where(abs(xlin - _x) < 0.0000001)
        yy = np.where(abs(ylin - _y) < 0.0000001)
        karr[xx, yy] = _z + 0.00000001

    condition = karr > 0.0000000001
    karrnonzero = np.extract(condition, karr)

    for _x, _y, _z in zip(x, y, z):
        xx = np.where(abs(xlin - _x) < 0.0000001)
        yy = np.where(abs(ylin - _y) < 0.0000001)
        karr2[xx, yy] = _z 

    return karr2[xi:xf, yi:yf], condition[xi:xf, yi:yf]


def calc_hist2d(A, nw0, nw1, condition):
    Nmax = A.shape
    N0 = int((Nmax[0] - (Nmax[0] % nw0)) / nw0)
    N1 = int((Nmax[1] - (Nmax[1] % nw1)) / nw1)
    k = np.zeros((N0, N1))
    kcond = np.zeros((N0, N1))
    for i in range(0, N0):
        ihead = (i+1)*nw0 - 1
        for j in range(0, N1):
            jhead = (j+1)*nw1 - 1
            if i == 0 and j == 0:
                k[i, j] = A[ihead, jhead]
            elif j == 0 and i != 0:
                k[i, j] = A[ihead, jhead] - A[ihead - nw0, jhead]
            elif i == 0 and j != 0:
                k[i, j] = A[ihead, jhead] - A[ihead, jhead - nw1]
            else:
                k[i, j] = A[ihead, jhead] - A[ihead - nw0, jhead] - A[ihead, jhead - nw1] + A[ihead - nw0, jhead - nw1]
    for i in range(0, N0):
        for j in range(0, N1):
            kcond[i, j] = np.sum(condition[i*nw0:(i+1)*nw0, j*nw1:(j+1)*nw1])
    return k, kcond


def calc_cost2d(A, maxw, condition):
    Cn = np.zeros((maxw))
    kaves = np.zeros((maxw))
    deltas = np.zeros((maxw))
    for i in range(1, maxw[0]):
        for j in range(1, maxw[1]):
            k, kcond = calc_hist2d(A, i, j, condition)
            knonzero = np.extract(np.max(kcond) == kcond, k)
            kave = np.average(knonzero)
            v = np.var(knonzero)
            cost = (2 * kave - v) / ((i*j)**2*1.0)
            Cn[i, j] = cost
            kaves[i, j] = kave
            deltas[i, j] = (i*j*1.0)
    return Cn, kaves, deltas


def get_optindx(m, n, Cn, kaves, deltas):
    ex = np.zeros((deltas.shape[0], deltas.shape[1]))
    ex[1:, 1:] = (1/m - 1/n) * kaves[1:, 1:] / (deltas[1:, 1:]**2*n) 
    ex[0, :] = 0.0
    ex[:, 0] = 0.0
    Cm = ex + Cn
    #print("opt bin index", np.unravel_index(np.argmin(Cm, axis=None), Cm.shape), "for m = ", m,  " with n = ", n)
    opt_indx = np.unravel_index(np.argmin(Cm, axis=None), Cm.shape)
    return opt_indx


def runex():
    #head = "/home/kazu/desktop/191031/"
    head = "/home/kazu/desktop/191120/"
    xi = 110
    xf = 217
    yi = 100
    yf = 220
    num_txtfiles = 16
    outfile = "result.txt"
    TotalIntensity = []
    opt_indx_for_x = []
    opt_indx_for_y = []
    TotalIntensity_ex = []
    for i in range(1, num_txtfiles + 1):
        txtfile = head + str(i) + "h.txt"
        #data, condition = get2ddata(txtfile, xi, xf, yi, yf)
        data, condition = get2ddata_for_commandline(txtfile, xi, xf, yi, yf)
        TotalIntensity.append(np.sum(data)*1.0)
    print("I have obtained info of total intensities")

    for Indx_of_txtfile in range(1, num_txtfiles + 1):
        print("Extraporation with n = ", TotalIntensity[Indx_of_txtfile-1])
        local_opt_indx_for_x = []
        local_opt_indx_for_y = []
        local_TotalIntensity_ex = []
        txtfile = head + str(Indx_of_txtfile) + "h.txt"
        #data, condition = get2ddata(txtfile, xi, xf, yi, yf)
        data, condition = get2ddata_for_commandline(txtfile, xi, xf, yi, yf)
        n = np.sum(data)*1.0
        maxxwidth = np.min(np.sum(condition, axis=0)) // 2
        maxywidth = np.min(np.sum(condition, axis=1)) // 2
        maxw = np.array([maxxwidth, maxywidth])
        cumdata = np.cumsum(np.cumsum(data, axis=0), axis=1)

        Cn, kaves, deltas = calc_cost2d(cumdata, maxw, condition)
        Cn = Cn / (n**2)   # This is according to the Cn in NeCo(2007)

        for m in TotalIntensity:
            opt_indx = get_optindx(m, n, Cn, kaves, deltas)
            local_opt_indx_for_x.append(opt_indx[0])
            local_opt_indx_for_y.append(opt_indx[1])
            local_TotalIntensity_ex.append(m)
            tmp_opt_x = opt_indx[0]
            tmp_opt_y = opt_indx[1]
        for mn_ratio in np.arange(1.1, 20.1, 0.1):
            m = mn_ratio * TotalIntensity[Indx_of_txtfile - 1]
            opt_indx = get_optindx(m, n, Cn, kaves, deltas)
            if opt_indx[0] < tmp_opt_x or opt_indx[1] < tmp_opt_y:
                tmp_opt_x = opt_indx[0]
                tmp_opt_y = opt_indx[1]
                local_opt_indx_for_x.append(opt_indx[0])
                local_opt_indx_for_y.append(opt_indx[1])
                local_TotalIntensity_ex.append(m)
        opt_indx_for_x.append(local_opt_indx_for_x)
        opt_indx_for_y.append(local_opt_indx_for_y)
        TotalIntensity_ex.append(local_TotalIntensity_ex)
    with open(outfile, mode='w') as f:
        f.write("optimization results for 373 data: n, opt_indx_x, opt_indx_y, 1/n, 1/(opt_idnx_x*opt_indx_y) \n")
        for Indx_of_txtfile in range(0, num_txtfiles): 
            f.write("%e %d %d %e %e\n" %
             (
              TotalIntensity[Indx_of_txtfile],
              opt_indx_for_x[Indx_of_txtfile][Indx_of_txtfile],
              opt_indx_for_y[Indx_of_txtfile][Indx_of_txtfile],
              1/(TotalIntensity[Indx_of_txtfile]*1.0),
              1/(opt_indx_for_x[Indx_of_txtfile][
              Indx_of_txtfile]*opt_indx_for_y[Indx_of_txtfile][
              Indx_of_txtfile]*1.0)
             )
            )
        f.write("extraporation results for 373 data: m, opt_indx_x, opt_indx_y, 1/m, 1/(opt_idnx_x*opt_indx_y) \n")
        for Indx_of_txtfile in range(0, num_txtfiles): 
            f.write("For n = %e \n" % TotalIntensity[Indx_of_txtfile])
            #for Indx_of_TotalIntensity in range(0, num_txtfiles):
            for Indx_of_TotalIntensity in range(0, len(TotalIntensity_ex[Indx_of_txtfile])):
                f.write("%e %d %d %e %e\n" %
                 (
                  TotalIntensity_ex[Indx_of_txtfile][Indx_of_TotalIntensity],
                  opt_indx_for_x[Indx_of_txtfile][Indx_of_TotalIntensity],
                  opt_indx_for_y[Indx_of_txtfile][Indx_of_TotalIntensity],
                  1/(TotalIntensity_ex[Indx_of_txtfile][Indx_of_TotalIntensity]*1.0),
                  1/(opt_indx_for_x[Indx_of_txtfile][
                     Indx_of_TotalIntensity]*opt_indx_for_y[Indx_of_txtfile][
                     Indx_of_TotalIntensity]*1.0)
                 )
                )

    
    fig=plt.figure(figsize=(12, 20))
    
    xlist_for_each_n = []
    ylist_for_each_n = []
    for Indx_of_txtfile in range(0, num_txtfiles):
        print(1/TotalIntensity[Indx_of_txtfile]*1.0)
        xlist_for_each_n.append(1/TotalIntensity[Indx_of_txtfile]*1.0)
        ylist_for_each_n.append(1/(opt_indx_for_x[Indx_of_txtfile][
               Indx_of_txtfile]*opt_indx_for_y[Indx_of_txtfile][
               Indx_of_txtfile]*1.0))
    for Indx_of_txtfile in range(0, num_txtfiles):
        xlist_for_m = []
        ylist_for_m = []
        for Indx_of_TotalIntensity in range(Indx_of_txtfile, len(TotalIntensity_ex[Indx_of_txtfile])):
            xlist_for_m.append(1/(TotalIntensity_ex[Indx_of_txtfile][Indx_of_TotalIntensity]*1.0))
            ylist_for_m.append(1/(opt_indx_for_x[Indx_of_txtfile][
                        Indx_of_TotalIntensity]*opt_indx_for_y[Indx_of_txtfile][
                        Indx_of_TotalIntensity]*1.0))
        if (Indx_of_txtfile+1 <= num_txtfiles/2):
            ax=fig.add_subplot(num_txtfiles//2, 2, 2*(Indx_of_txtfile+1)-1)
        else:
            ax=fig.add_subplot(num_txtfiles//2, 2, 2*(Indx_of_txtfile+1)-num_txtfiles)
        ax.scatter(xlist_for_each_n[Indx_of_txtfile:],
                ylist_for_each_n[Indx_of_txtfile:], marker='x', clip_on=False,
                s=50, label="each n")
        ax.scatter(xlist_for_m, ylist_for_m,
                marker='+', clip_on=False, s=72, label="prediction")
        ax.set_xlim(0,0.0001)
        ax.set_ylim(0,0.06)
        hour = Indx_of_txtfile + 1 
        ax.text(0.25, 0.8, 'n at %dh'%hour,
        transform=ax.transAxes, ha="right")

        if (Indx_of_txtfile+1 == 1): 
           ax.legend()

        if (Indx_of_txtfile+1 == 1 or Indx_of_txtfile+1 == num_txtfiles//2 + 1): 
           ax.set_yticks([0,0.02,0.04,0.06])
        else:
           ax.set_yticks([0,0.02,0.04])

        ax.tick_params(labelbottom=False)
        ax.tick_params(direction = "in")
        if (Indx_of_txtfile+1 ==  num_txtfiles//2  or Indx_of_txtfile+1 == num_txtfiles):
            ax.tick_params(labelbottom=True)
            ax.set_xlabel('1/m or 1/n')

        #plt.gca().ticklabel_format(style="sci", scilimits=(0,0), axis="y")
        plt.gca().ticklabel_format(style="sci", scilimits=(0,0), axis="x")
        if (Indx_of_txtfile == num_txtfiles//4 or Indx_of_txtfile == num_txtfiles//2 + num_txtfiles//4 ):
            ax.set_ylabel('1/(opt_wx*opt_wy)')
    plt.subplots_adjust(wspace=0.4, hspace=0.0)
    #plt.savefig("result.txt.full.pdf")
    plt.show()

      

runex()


