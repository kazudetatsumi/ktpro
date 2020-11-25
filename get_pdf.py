#!/usr/bin/env python
# This script generates a probability density function and store it as a h5py file.
# Prior to execute this script, double differential cross-sections should be calculated and stored in a h5py file by using a script calc_ddscs_BL*.py.
# The pdf is generated by convolving the calculated ddscs with appropriate point spread functions, adjusted to the beam-lines. 
# Kazuyoshi TATSUMI 2019.
import numpy as np
import h5py
import yaml
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from scipy import *
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter

fig = plt.figure(figsize=(9, 9))
plt.rcParams['font.family'] = 'Arial'
fig.suptitle("crosssections of 4D INS phantom data for Ei42", fontsize="x-large")

Ei = 50.  # Incident Neutron Energy [meV]

def store_omega_in_goodarray(fn, mesh):
    f = h5py.File(fn, 'r')
    nmodes = f["frequency"].shape[1]      # f["frequency"]:[number of qpoints, number of modes], ordered as z indx increasing faster
    omega = np.reshape(f["frequency"][:],
                       [mesh[0], mesh[1], mesh[2], nmodes], order='C')  
    omega[omega < 0] = 0.0
    print("omega readed from qpoints.hdf5")
    return omega


def get_pdf(omega, gamma, de,  ddscs):
    #print("max energy", np.max(omega))
    #ene = np.arange(np.min(omega), np.max(omega)+de, de)
    ene = np.linspace(0.0, 40.0, int((40.0 + 0.0) // de) + 1)
    print(ene[0:5], "<- check energy range")
    print(ene[ene.shape[0]-5:], "<- check energy range")
    nx = omega.shape[0]
    ny = omega.shape[1]
    nz = omega.shape[2]
    no = omega.shape[3]
    nene = ene.shape[0]
    pdf = np.zeros((nx, ny, nz, nene))  # pdf = probability density function
    for i in range(0, nx):
        for j in range(0, ny):
            for k in range(0, nz):
                for h in range(0, no):
                    if omega[i, j, k, h] < Ei:
                        #pdf[i, j, k, :] += fun_lore(ene, gamma, omega[i, j, k, h]) * ddscs[i, j, k, h]
                        pdf[i, j, k, :] += fun_gaus(ene, omega[i, j, k, h]) * ddscs[i, j, k, h]
    for eidx in range(0, pdf.shape[3]):
        pdf[:, :, :, eidx] = gaussian_filter(pdf[:, :, :, eidx], sigma=[2.3/2.3548, 2.3/2.3548, 2.3/2.3548])

    pdf = pdf / ((np.sum(pdf)))
    print("pdf is created.")
    return pdf


def fun_lore(x, gamma, mean):
    dx = x[1] - x[0]
    l = (1 / 3.14) * 0.5 * gamma / ((x - mean)**2 + (0.5 * gamma)**2)
    l = l/(np.sum(l)*dx)
    return l


def fun_gaus(x, mean):
    sigma = eres_siki(mean) / 2.3548
    dx = x[1] - x[0]
    l = 1.0 / (sigma * (2.0*np.pi)**0.5) * np.exp(-0.5 *
                                                  ((x - mean) / sigma)**2)
    l = l/(np.sum(l)*dx)
    return l


def get_dtm(Ei):
    Eis = [
        0.099999, 0.12589, 0.15849, 0.19953, 0.25119, 0.31623, 0.3981, 0.50118,
        0.63096,  0.79433, 0.99999, 1.2589,  1.5849,  1.9953,  2.5119, 3.1623,
        3.981,    5.0118,  6.3096,  7.9433,  9.9999,  12.589,  15.849, 19.953,
        25.119,   31.623,  39.81,   50.118,  63.096,  79.433,  99.999, 125.89,
        158.49,   199.53,  251.19,  316.23,  398.1,   501.18,  630.96, 794.33,
        999.99,   1258.9,  1584.9,  1995.3,  2511.9,  3162.3,  3981,   5011.8,
        6309.6,   7943.3,  9999.9
        ]
    Dtms = [
        487.05,   487.71,  461.9,   448,     417.73,  383.09,  363.51, 337.71,
        307.18,   286.35,  262.2,   241.77,  227.26,  208.89,  188.66, 172.44,
        155.93,   137.95,  119.63,  105.82,  89.559,  76.105,  58.736, 42.069,
        31.527,   25.151,  20.432,  17.571,  15.06,   13.605,  12.149, 10.844,
        9.6387,   8.5343,  7.5804,  6.7772,  5.974,   5.472,   4.8193, 4.4679,
        4.0161,   3.5141,  3.1627,  2.8113,  2.5603,  2.2591,  2.0583, 1.8073,
        1.6566,   1.4558,  1.3052
        ]
    Dtm0 = interp1d(Eis, Dtms)
    return(Dtm0(Ei))


def eres_ama(mean):
    L1 = 30.0  # [m]
    L2 = 4.0   # [m]
    L3 = 1.6   # [m]
    # Ei = 42.05  # [meV]
    # Ei = 23.65 # [meV]
    tch = (5227000.0/Ei)**0.5*(L1-L3)
    wslit = 30.0  # [mm]
    ddisk = 600.0  # [mm]
    f = 150.0   # [Hz]
    # f = 300.0   # [Hz]
    dtch = wslit/(ddisk*np.pi*f)*10.0**6/2.0  # [micro_sec]
    # dtm  = 2.0/(Ei*0.001)**0.5 # [micro_sec]
    dtm = get_dtm(Ei)
    # ssize = 20.0 # [mm]
    # wdet =  0.019 # [m]
    # dL2 = (wdet**2 + ssize**2)**0.5*0.001 # [m]
    # dL2 = sqrt((wdet/4.0*3.1415927)**2+(ssize/1000./4.0*3.1415927)**2)

    first = dtch/tch*(1.0+L1/L2*(1.0-mean/Ei)**1.5)
    second = dtm/tch*(1.0+L3/L2*(1.0-mean/Ei)**1.5)
    # third = dL2/L2*(1.0-mean/Ei)
    deofaei = 2.*(first**2. + second**2.)**0.5
    return(deofaei*Ei)


def eres_siki(mean):
    L1 = 18.0  # moderator-sample distance [m]
    L2 = 2.5   # sample-detector distance [m]
    L3 = 1.71  # chopper-sample distance [m]
    wdet = 0.019  # width of a detector tube [m]

    Lgs = 2200.     # guide-sample distance [mm]
    wguide = 43.35  # width of the end of the guide [mm]
    mvalue = 4.     # m-value of the supermirror
    gslope = (56.74-43.35)/2./1755.*1000        # slope of the last section of the guide [mrad]
    wslit = 0.44    # width of each slit of the chopper [mm]
    lslit = 20.0  # length of each slit of the chopper [mm]
    div_chopper = wslit/lslit*1000.   # angular divergence defined by the chopper [mrad]

    ssize = 20  # Sample size [mm]
    freq = 150  # Chopper frequency [Hz]
    #Ei = 50.  # Incident Neutron Energy [meV]

    # pulse width at the moderator [usec]
    Dtm = get_dtm(Ei)

    # time when neutrons reach the chopper [usec]
    tch = 2286. / sqrt(Ei) * (L1-L3)

    # angular divergence of the incident beam [mrad]
    #  (1) divergence by reflection by the supermirror guide:
    #        slope of the last section of the guide
    #     + critical angle of the mirror
    div_guide = 2.*(gslope + 1.73*mvalue*9.044/sqrt(Ei))
    #  (2) divergence defined by the sample [mrad]
    div_sample = 2.*(wguide/2. + ssize/2.) / Lgs * 1000.
    div_beam = min(div_guide, div_sample)

    # broadening of the pulse width at the chopper
    # due to the light-house effect [cf. Windsor's textbook]
    u = div_beam/div_chopper
    if u <= 0.8:
        p = 1. + u/4.
    else:
        if u <= 2.:
            p = 2. + u - sqrt(4.*u-u*u)
        else:
            p = u

    # opening time of the chopper [usec]
    Dtch0 = wslit/(2.*pi*lslit*freq)*1e6   # nominal value
    Dtch = Dtch0 * p   # broadening due to the lighthouse effect

    # uncertinity of L2 due to the sample size and the detector size
    DL2 = sqrt((wdet/4.0*3.1415927)**2+(ssize/1000./4.0*3.1415927)**2)

    # energy resolution relative to Ei [cf. Windsor's textbook]
    Dmean_Ei = sqrt(
                (2. * Dtch/tch * (1. + L1/L2 * (1.-mean/Ei)**1.5))**2 +
                (2. * Dtm/tch * (1. + L3/L2 * (1.-mean/Ei)**1.5))**2 +
                (2. * DL2/L2 * (1. - mean/Ei))**2
               )
    return(Dmean_Ei*Ei)


def save_h5py(data4, outfile):
    with h5py.File(outfile, 'w') as hf:
        hf.create_dataset('data4', data=data4)
       

def read_h5py(outfile, term):
    f=h5py.File(outfile, 'r')
    return np.array(f[term])


def save_pdf_hdf5(elimifile, pdf, condition):
    with h5py.File(elimifile, 'w') as hf:
        hf.create_dataset('pdf', data=pdf)
        hf.create_dataset('condition', data=condition)


def plotter(xi, xe, yi, ye, lx, ly, data, vn, hn, cn, dev, text):
    ax = fig.add_subplot(vn,  hn, cn)
    ax.pcolor(np.transpose(data), vmax=np.max(data[xi:xe, yi:ye])/dev, cmap='jet')
    ax.set_xlabel(lx)
    ax.set_ylabel(ly)
    ax.axvline(x=xi, color='white', lw=0.5)
    ax.axvline(x=xe, color='white', lw=0.5)
    ax.axhline(y=yi, color='white', lw=0.5)
    ax.axhline(y=ye, color='white', lw=0.5)
    ax.xaxis.set_label_coords(0.5, 1.145)
    ax.tick_params(direction="in", color="white", top=True, labeltop=True, labelbottom=False)
    ax.text(2, 75, text, color='white')
    ax.axis('tight')


def plot_crosssection_more_useful(xi, xe, yi, ye, zi, ze, ei, ee, condition, data4):
    print(np.sum(data4))
    plotter(xi, xe, yi, ye, 'qx', 'qy', np.sum(data4[:, :, :, 0], axis=2),     3, 2, 1, 1.0, "qz:all, E=0")
    plotter(xi, xe, yi, ye, 'qx', 'qy', np.sum(condition[:, :, :, 0], axis=2), 3, 2, 2, 1.0, "")
    plotter(yi, ye, ei, ee, 'qy', 'E',  data4[106, :, 34, :],                   3, 2, 3, 1.0, "qx=106, qz=34")
    plotter(yi, ye, ei, ee, 'qy', 'E',  condition[106, :, 34, :],               3, 2, 4, 1.0, "")
    plotter(xi, xe, ei, ee, 'qx', 'E',  data4[:, 146, 34, :],                    3, 2, 5, 1.0, "qy=146, qz=34")
    plotter(xi, xe, ei, ee, 'qx', 'E',  condition[:, 146, 34, :],                3, 2, 6, 1.0, "")


def run():
    dqx = 0.025
    dqy = 0.025
    dqz = 0.025
    lb_qx = -1.65
    ub_qx = 4.10
    lb_qy = -2.1
    ub_qy = 2.8
    lb_qz = -0.85
    ub_qz = 0.9
    dq = np.array([dqx, dqy, dqz])
    lb_q = np.array([lb_qx, lb_qy, lb_qz])
    ub_q = np.array([ub_qx, ub_qy, ub_qz])
    xlin = np.arange(lb_qx, ub_qx, dqx)
    nx = xlin.shape[0]
    ylin = np.arange(lb_qy, ub_qy, dqy)
    ny = ylin.shape[0]
    zlin = np.arange(lb_qz, ub_qz, dqz)
    nz = zlin.shape[0]
    mesh = np.array([nx, ny, nz])
    print(mesh)
    print(mesh[0]*mesh[1]*mesh[2])

    PlanckConstant = 4.13566733e-15       # [eV s]
    THzTomev = PlanckConstant * 1e15      # [meV]
    gamma = 0.8                           # FWHM [meV]
    dE = 0.5                              # [meV]
    qpointsfile = "/home/kazu/WORK/vasp-phonopy/cu_BL01/test/qpoints.hdf5"
    ddscsfile = "/home/kazu/WORK/vasp-phonopy/cu_BL01/test/ddscs.hdf5"
    maskfile = "/home/kazu/desktop/200204/fine/out_hw_all.hdf5"
    ## For first run, uncomment the three lines below.
    #omega = store_omega_in_goodarray(qpointsfile, mesh)*THzTomev  # [meV]
    #ddscs = read_h5py(ddscsfile, "ddscs")
    #pdf = get_pdf(omega, gamma, dE, ddscs)
    ## For first run, comment out the two lines below.
    pdf = read_h5py("pdf.hdf5", "pdf")
    condition = read_h5py("pdf.hdf5", "condition")
    xi = 120
    xe = 172
    yi = 61
    ye = 145
    zi = 16
    ze = 53
    ei = 20
    ee = 70
    #condition = read_h5p(maskfile, condition)
    #save_pdf_hdf5("pdf.hdf5", pdf[:, :, :, :], condition[:, :, :, :])
    plot_crosssection_more_useful(xi, xe, yi, ye, zi, ze, ei, ee, condition, pdf)
    plt.savefig("pdf.png")


run()
#plt.show()
