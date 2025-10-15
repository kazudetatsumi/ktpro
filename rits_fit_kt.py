import sys
import os
import numpy as np
import scipy.special as sp
import scipy.optimize as so

np.random.seed(222)

def get_sim_spectrum(inpfile='rits_initial.inp'):
    sys.path.append("/home/kazu/desktop/240424/rits/src/cpp/rietveldcc/Linux")
    from AdvRits import AdvRietveld
    rits = AdvRietveld()
    if not os.path.exists("spgra"):
        print("'spgra' dose not exist.")
        return

    rits.SetInitialFileName(inpfile)
    rits.ReadInitialFile()
    rits.SetDataFileName('fit_results.dat.tmp')
    rits.ReadDataFile()
    rits.SetFittingResultDataFileName('fit_func.data')
    rits.SetFittingResultParamFileName('fit_para.data')

    rits.Fit()

    return np.array(rits.PutObservedLambda()),\
        np.array(rits.PutFittedIntensity())


def get_sim_edgespectrum(inpfile='edge_3.inp'):
    sys.path.append("/home/kazu/desktop/240424/rits/src/cpp/edgecc/Linux")
    from AdvEdge import AdvEdge
    rits = AdvEdge()
    rits.SetInitialFileName(inpfile)
    rits.ReadInitialFile()
    rits.SetDataFileName('fit_results.dat.tmp')
    rits.ReadDataFile()
    rits.SetFittingResultDataFileName('fit_func.data')
    rits.SetFittingResultParamFileName('fit_para.data')

    rits.Fit()

    return np.array(rits.PutObservedLambda()),\
        np.array(rits.PutFittedIntensity())


def get_sim_edgespectrum_local2(para):
    flight = para[15]
    thkl = 2.*para[4]*flight*1.0E6/3956.
    sigm, alph, beta = setFacilityParameter(para)
    tofIn = np.arange(152)*20+23000
    delt = tofIn - thkl
    Pu = 0.5 * alph * (alph*sigm*sigm + 2.*delt)
    Pv = 0.5 * beta * (beta*sigm*sigm - 2.*delt)
    Py = (alph*sigm*sigm + delt) / (np.sqrt(2.) * sigm)
    Pz = (beta*sigm*sigm - delt) / (np.sqrt(2.) * sigm)
    Pw = delt / (np.sqrt(2.) * sigm)
    # I heared that the latest rits uses double precision and approximationexperfc is not needed.
    #exp_erf1 = ApproximationExpErfc(Py, Pu)
    #exp_erf2 = ApproximationExpErfc(Pz, Pv)
    #step = 0.5*sp.erfc(Pw) - 0.5 * (beta*exp_erf1 - alph * exp_erf2) / (alph + beta)
    step = 0.5*sp.erfc(Pw) - 0.5 * (beta*np.exp(Pu)*sp.erfc(Py) - alph * np.exp(Pv)*sp.erfc(Pz)) / (alph + beta)
    yCalc = np.exp((-1.) * ((step*(para[2] + para[3] * (1.0E-4)*tofIn)) + (para[0] + para[1] * (1.0E-4)*tofIn)))
    return tofIn, yCalc


def get_sim_edgespectrum_local(inpfile='edge_3.inp'):
    para = get_para(inpfile)
    return get_sim_edgespectrum_local2(para)


def err_transmission(sample, openbeam, fac):
    err_sample = sample**0.5
    err_openbeam = openbeam**0.5
    return ((1./openbeam*err_sample)**2 +
            (sample/openbeam**2*err_openbeam)**2)**0.5*fac


def test_fit(inpfile):
    import pickle
    with open('bi3d_scratch_255_partial_phantom_local_mean.pkl', 'rb') as f:
        sample = pickle.load(f).transpose((0, 2, 1))
        op = pickle.load(f).transpose((0, 2, 1))
    with open('../params_scratch.pkl', 'rb') as f:
        params = pickle.load(f)[:, ::-1, :]
    for i in range(4):
        params[i][params[i]!=0] = np.round(np.mean(params[i][params[i]!=0]), 2)
    para = get_para(inpfile)
    variables = para[4:6]
    #print(sample.shape)
    for vpos in range(sample.shape[1]):
        #hpos=100
        hpos=31
        y = ((sample/op)[:, vpos, hpos], err_transmission(sample[:, vpos, hpos], op[:, vpos, hpos], 1.))
        para[:4] = params[:4, vpos, hpos]
    #checkres(variables, para, y)
    #results = so.least_squares(res, variables, loss='linear', f_scale=0.1, args=(para, y))
        results = so.least_squares(res, variables, args=(para, y))
        print('vpos=', vpos, 'optimized')
        #print(results.message)
        if vpos == 0:
            results_x = results.x
        else:
            results_x = np.vstack((results_x, results.x))
    print(results_x.shape)
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(3, 1)
    ax[0].plot(np.abs(results_x[:, 1]))
    ax[0].plot(params[5, :, hpos])
    #plt.plot(np.sum((sample/op)[:, vpos, :], axis=0))
    #ax[0].set_ylim([np.min(results_x[:, 0]), np.max(params[4, vpos])])
    ax[1].imshow(np.sum(sample/op, axis=0), origin='lower')
    #ax[1].imshow(params[4], origin='lower')
    ax[2].imshow(params[4, :, :], origin='lower')
    plt.show()
    #print(results)
    #checkres(results.x, para, y)
    #Jacobian = results.jac
    #cov_matrix = np.linalg.inv(Jacobian.T @ Jacobian)
    #errorbars = np.sqrt(np.diag(cov_matrix))
    #print(errorbars)


def checkres(variables, para, y):
    import matplotlib.pyplot as plt
    para[4] = variables[0]
    para[5] = variables[1]
    tofIn, yCalc = get_sim_edgespectrum_local2(para)
    plt.plot(yCalc)
    plt.plot(y[0])
    plt.ylim([0, 1])
    plt.show()


def res(variables, para, y):
    para[4] = variables[0]
    para[5] = variables[1]
    tofIn, yCalc = get_sim_edgespectrum_local2(para)
    return (yCalc - y[0])/y[1]


def get_para(inpfile):
    para = []
    for il, line in enumerate(open(inpfile)):
        if "?" in line:
            para.append(float(line.split()[0][1:]))
        else:
            para.append(float(line.split()[0]))
    return np.array(para)


def ApproximationExpErfc_(x, y):
    CfOrder = 20
    CfThres = 2.4
    func = np.zeros_like(x) + 0.01
    area1 = np.abs(x) < CfThres
    area2 = (np.abs(x) > 1.0e50) & (x > 0.)
    area3 = (np.abs(x) > 1.0e50) & (x <= 0.)
    area4 = (np.abs(x) <= 1.0e50) & (np.abs(x) >= CfThres) & (x >= 0.)
    area5 = (np.abs(x) <= 1.0e50) & (np.abs(x) >= CfThres) & (x < 0.)
    func[area1] = np.exp(y[area1])*sp.erfc(x[area1])
    func[area2] = 0.
    func[area3] = -2.*np.exp(y[area3])
    d = np.abs(x[area4])*np.sqrt(2.)
    a = np.zeros_like(x[area4])
    for n in np.arange(CfOrder, 0, -1):
        a = n / (d + a)
    a = np.sqrt(2./np.pi)*np.exp(y[area4]-x[area4]*x[area4])/(d+a)
    func[area4] = a
    d = np.abs(x[area5])*np.sqrt(2.)
    a = np.zeros_like(x[area5])
    for n in np.arange(CfOrder, 0, -1):
        a = n / (d + a)
    a = np.sqrt(2./np.pi)*np.exp(y[area5]-x[area5]*x[area5])/(d+a)
    func[area5] = 2.*np.exp(y[area5]) - a
    return func


def ApproximationExpErfc(x, y):
    CfOrder = 20
    CfThres = 2.4
    func = np.zeros_like(x) + 0.01
    area1 = np.abs(x) < CfThres
    area2 = (np.abs(x) > 1.0e50) & (x > 0.)
    area3 = (np.abs(x) > 1.0e50) & (x <= 0.)
    area4 = (np.abs(x) <= 1.0e50) & (np.abs(x) >= CfThres)
    area41 = (np.abs(x) <= 1.0e50) & (np.abs(x) >= CfThres) & (x >= 0.)
    area42 = (np.abs(x) <= 1.0e50) & (np.abs(x) >= CfThres) & (x < 0.)
    _area41 = x[area4] >= 0.
    _area42 = x[area4] < 0.
    func[area1] = np.exp(y[area1])*sp.erfc(x[area1])
    func[area2] = 0.
    func[area3] = -2.*np.exp(y[area3])
    d = np.abs(x[area4])*np.sqrt(2.)
    a = np.zeros_like(x[area4])
    for n in np.arange(CfOrder, 0, -1):
        a = n / (d + a)
    a = np.sqrt(2./np.pi)*np.exp(y[area4]-x[area4]*x[area4])/(d+a)
    func[area41] = a[_area41]
    func[area42] = 2.*np.exp(y[area42]) - a[_area42]
    return func


def setFacilityParameter(para):
    sigm1 = para[5] * 10.
    sigm2 = para[6] * 10.
    alph = para[7] * 0.01
    beta = para[8] * 0.01
    sigm = np.sqrt(sigm1**2 + sigm2**2)
    return sigm, alph, beta





#x, y = get_sim_spectrum()
#x = np.array(x)
#y = np.array(y)

#timescale = 3600000
#y *= timescale
#sensitivity = 0.01
#import matplotlib.pyplot as plt
#_y = np.random.poisson(y*sensitivity*x)/(sensitivity*x)
#plt.scatter(x, _y, s=15, fc='none', ec='b')
#plt.plot(x, y, c='r')
#plt.plot(x, y - _y) 
#plt.xlim([1, 5.1])
#plt.show()

test_fit('edge_3.inp.phantom')
