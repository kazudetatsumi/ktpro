import numpy as np
import scipy.special as sp
import scipy.optimize as so


def get_tofIn():
    return np.arange(152)*20+23000


def get_thkl(para):
    flight = para[15]
    return 2.*para[4]*flight*1.0E6/3956.


def get_postedgearea(para):
    tofIn = get_tofIn()
    thkl = get_thkl(para)
    return (tofIn >= para[13]*thkl) & (tofIn <= para[14]*thkl)


def get_preedgearea(para):
    tofIn = get_tofIn()
    thkl = get_thkl(para)
    return (tofIn >= para[9]*thkl) & (tofIn <= para[10]*thkl)


def get_edgearea(para):
    tofIn = get_tofIn()
    thkl = get_thkl(para)
    return (tofIn >= para[11]*thkl) & (tofIn <= para[12]*thkl)


def get_fulledgearea(para):
    tofIn = get_tofIn()
    thkl = get_thkl(para)
    return (tofIn >= para[9]*thkl) & (tofIn <= para[14]*thkl)


def get_allarea(para):
    tofIn = get_tofIn()
    return (tofIn >= np.min(tofIn)) & (tofIn <= np.max(tofIn))


def get_sim_edgespectrum_local(para):
    thkl = get_thkl(para)
    sigm, alph, beta = setFacilityParameter(para)
    tofIn = get_tofIn()
    delt = tofIn - thkl
    Pu = 0.5 * alph * (alph*sigm*sigm + 2.*delt)
    Pv = 0.5 * beta * (beta*sigm*sigm - 2.*delt)
    Py = (alph*sigm*sigm + delt) / (np.sqrt(2.) * sigm)
    Pz = (beta*sigm*sigm - delt) / (np.sqrt(2.) * sigm)
    Pw = delt / (np.sqrt(2.) * sigm)
    step = 0.5*sp.erfc(Pw) - 0.5 * (beta*np.exp(Pu)*sp.erfc(Py) - alph *
                                    np.exp(Pv)*sp.erfc(Pz)) / (alph + beta)
    yCalc = np.exp((-1.) * ((step*(para[2] + para[3] * (1.0E-4)*tofIn)) +
                            (para[0] + para[1] * (1.0E-4)*tofIn)))
    return yCalc


def get_sim_edgespectrum_local3(para, edgearea):
    thkl = get_thkl(para)
    tofIn = get_tofIn()[edgearea]
    sigm, alph, beta = setFacilityParameter(para)
    delt = tofIn - thkl
    Pu = 0.5 * alph * (alph*sigm*sigm + 2.*delt)
    Pv = 0.5 * beta * (beta*sigm*sigm - 2.*delt)
    Py = (alph*sigm*sigm + delt) / (np.sqrt(2.) * sigm)
    Pz = (beta*sigm*sigm - delt) / (np.sqrt(2.) * sigm)
    Pw = delt / (np.sqrt(2.) * sigm)
    step = 0.5*sp.erfc(Pw) - 0.5 * (beta*np.exp(Pu)*sp.erfc(Py) - alph *
                                    np.exp(Pv)*sp.erfc(Pz)) / (alph + beta)
    yCalc = np.exp((-1.) * ((step*(para[2] + para[3] * (1.0E-4)*tofIn)) +
                            (para[0] + para[1] * (1.0E-4)*tofIn)))
    return yCalc


def get_sim_edgespectrum_local1(para, postedgearea):
    tofIn = get_tofIn()[postedgearea]
    return np.exp((-1.) * (para[0] + para[1] * (1.0E-4)*tofIn))


def get_sim_edgespectrum_local2(para, preedgearea):
    tofIn = get_tofIn()[preedgearea]
    return np.exp((-1.)*(para[0] + para[2] +
                  (para[1] + para[3])*(1.0E-4)*tofIn))


def err_transmission(sample, openbeam, fac):
    err_sample = sample**0.5
    err_openbeam = openbeam**0.5
    return ((1./openbeam*err_sample)**2 +
            (sample/openbeam**2*err_openbeam)**2)**0.5*fac


def get_bi3d_and_their_params():
    import pickle
    with open('bi3d_scratch_255_partial_phantom_local_mean.pkl', 'rb') as f:
        sample = pickle.load(f).transpose((0, 2, 1))
        op = pickle.load(f).transpose((0, 2, 1))
    with open('params_scratch_for_class.pkl', 'rb') as f:
        params = pickle.load(f)
    for i in range(4):
        params[i][params[i] != 0] = np.round(np.mean(params[i][params[i] != 0]
                                                     ), 2)
    y = (sample/op, err_transmission(sample, op, 1.))
    return y, params


def plot_results4(results_x, params, y, hpos):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(8, 1)
    for aidx, pidx in enumerate([4, 5, 0, 1, 2, 3]):
        ax[aidx].plot(results_x[:, pidx])
        ax[aidx].plot(params[pidx, :, hpos])
        ax[aidx].set_ylim([np.min(np.unique(results_x[:, pidx])[1:]),
                           np.max(np.unique(results_x[:, pidx])[1:])])
        ax[aidx].set_xlim([0, results_x.shape[0]])
    ax[6].imshow(np.sum(y[0], axis=0), origin='lower')
    ax[7].imshow(params[4, :, :], origin='lower')
    plt.show()


def plot_results3(results_x, params, y, hpos):
    import matplotlib.pyplot as plt
    print(results_x.shape)
    fig, ax = plt.subplots(4, 1)
    for aidx, pidx in enumerate([4, 5]):
        ax[aidx].plot(results_x[:, pidx-4])
        ax[aidx].plot(params[pidx, :, hpos])
        ax[aidx].set_ylim([np.min(np.unique(results_x[:, pidx-4])[1:]),
                           np.max(np.unique(results_x[:, pidx-4])[1:])])
        ax[aidx].set_xlim([0, results_x.shape[0]])
    ax[2].imshow(np.sum(y[0], axis=0), origin='lower')
    ax[3].imshow(params[4, :, :], origin='lower')
    plt.show()


def plot_results1(results_x, params, y, hpos):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(3, 1)
    ax[0].plot(results_x[:, 0])
    ax[0].plot(params[0, :, hpos])
    ax[1].imshow(np.sum(y[0], axis=0), origin='lower')
    ax[2].imshow(params[0, :, :], origin='lower')
    plt.show()


def plot_results12(results_x, params, y, hpos):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(3, 1)
    ax[0].plot(results_x[:, 1])
    ax[0].plot(params[3, :, hpos])
    ax[1].imshow(np.sum(y[0], axis=0), origin='lower')
    ax[2].imshow(params[3, :, :], origin='lower')
    plt.show()


def test_fit1(inpfile):
    y, params = get_bi3d_and_their_params()
    para = get_para(inpfile)
    variables = para[0:2]
    postedgearea = get_postedgearea(para)
    hpos = 31
    for vpos in range(y[0].shape[1]):
        results = so.least_squares(res1, variables,
                                   args=(para,
                                         (y[0][postedgearea, vpos, hpos],
                                          y[1][postedgearea, vpos, hpos]
                                          ), postedgearea))
        if vpos == 0:
            results_x = results.x
        else:
            results_x = np.vstack((results_x, results.x))
    plot_results1(results_x, params, y, hpos)


def test_fit12(inpfile):
    y, params = get_bi3d_and_their_params()
    para = get_para(inpfile)
    variables = para[0:2]
    postedgearea = get_postedgearea(para)
    preedgearea = get_preedgearea(para)
    hpos = 31
    for vpos in range(y[0].shape[1]):
        results = so.least_squares(res1, variables,
                                   args=(para, (y[0][postedgearea, vpos, hpos],
                                                y[1][postedgearea, vpos, hpos]
                                                ), postedgearea))
        para[0:2] = results.x
        variables = para[2:4]
        results = so.least_squares(res2, variables,
                                   args=(para, (y[0][preedgearea, vpos, hpos],
                                                y[1][preedgearea, vpos, hpos]
                                                ), preedgearea))
        para[2:4] = results.x
        if vpos == 0:
            results_x = results.x
        else:
            results_x = np.vstack((results_x, results.x))
    plot_results12(results_x, params, y, hpos)


def test_fit123(inpfile, numoftrials):
    y, params = get_bi3d_and_their_params()
    mask = get_mask(y[0])
    hpos = 30
    results_x = np.zeros((y[0].shape[1], 6))
    for vpos in range(y[0].shape[1]):
        if not mask[vpos, hpos]:
            print(vpos)
            para = get_para(inpfile)
            for nt in range(numoftrials):
                para = fit123(para, (y[0][:, vpos, hpos], y[1][:, vpos, hpos]))
            results_x[vpos] = para[:6]
    plot_results4(results_x, params, y, hpos)


def fit123(para, y):
    postedgearea = get_postedgearea(para)
    preedgearea = get_preedgearea(para)
    edgearea = get_edgearea(para)
    variables = para[0:2]
    results = so.least_squares(res1, variables, method='trf',
                               args=(para, (y[0][postedgearea],
                                            y[1][postedgearea]), postedgearea))
    para[0:2] = results.x
    variables = para[2:4]
    results = so.least_squares(res2, variables, method='trf',
                               args=(para, (y[0][preedgearea],
                                            y[1][preedgearea]), preedgearea))
    para[2:4] = results.x
    variables = para[4:6]
    results = so.least_squares(res3, variables, method='trf',
                               args=(para, (y[0][edgearea],
                                            y[1][edgearea]), edgearea))
    para[4:6] = results.x
    return para


def fit4(para, y):
    fulledgearea = get_fulledgearea(para)
    variables = para[:6]
    results = so.least_squares(res4, variables, method='trf',
                               args=(para, (y[0][fulledgearea],
                                            y[1][fulledgearea]), fulledgearea))
    para[:6] = results.x
    return para


def test_fit1234(inpfile, numoftrials):
    y, params = get_bi3d_and_their_params()
    mask = get_mask(y[0])
    hpos = 30
    results_x = np.zeros((y[0].shape[1], 6))
    for vpos in range(y[0].shape[1]):
        if not mask[vpos, hpos]:
            print(vpos)
            para = get_para(inpfile)
            for nt in range(numoftrials):
                para = fit123(para, (y[0][:, vpos, hpos], y[1][:, vpos, hpos]))
            para = fit4(para, (y[0][:, vpos, hpos], y[1][:, vpos, hpos]))
            results_x[vpos] = para[:6]
    plot_results4(results_x, params, y, hpos)


def get_mask(transmission):
    _transmission = np.average(transmission[90:], axis=0)
    return _transmission >= 0.95


def test_fit(inpfile):
    # Simple case where only the two parames[4:6] should be fitted
    # over the whole tof area.
    # The other fitting parameters (parames[0:4] is fixed as the
    # true values. 
    y, params = get_bi3d_and_their_params()
    mask = get_mask(y[0])
    para = get_para(inpfile)
    variables = para[4:6]
    hpos = 31
    results_x = np.zeros((y[0].shape[1], 2))
    for vpos in range(y[0].shape[1]):
        if not mask[vpos, hpos]:
            para[:4] = params[:4, vpos, hpos]
            results = so.least_squares(res, variables,
                                       args=(para, (y[0][:, vpos, hpos],
                                                    y[1][:, vpos, hpos])))
            results_x[vpos] = results.x
    plot_results3(results_x, params, y, hpos)
    ## for fitting error estimation
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
    tofIn, yCalc = get_sim_edgespectrum_local(para)
    plt.plot(yCalc)
    plt.plot(y[0])
    plt.ylim([0, 1])
    plt.show()


def checkres1(variables, para, y, postedgearea):
    import matplotlib.pyplot as plt
    para[0] = variables[0]
    para[1] = variables[1]
    print(variables)
    yCalc = get_sim_edgespectrum_local1(para, postedgearea)
    plt.plot(yCalc)
    plt.plot(y[postedgearea])
    #plt.ylim([0, 1])
    plt.show()


def res(variables, para, y):
    para[4] = variables[0]
    para[5] = variables[1]
    yCalc = get_sim_edgespectrum_local(para)
    return (yCalc - y[0])/y[1]


def res4(variables, para, y, fulledgearea):
    para[:6] = variables
    yCalc = get_sim_edgespectrum_local3(para, fulledgearea)
    return (yCalc - y[0])/y[1]


def res3(variables, para, y, edgearea):
    para[4] = variables[0]
    para[5] = variables[1]
    yCalc = get_sim_edgespectrum_local3(para, edgearea)
    return (yCalc - y[0])/y[1]


def res2(variables, para, y, preedgearea):
    para[2] = variables[0]
    para[3] = variables[1]
    yCalc = get_sim_edgespectrum_local2(para, preedgearea)
    return (yCalc - y[0])/y[1]


def res1(variables, para, y, postedgearea):
    para[0] = variables[0]
    para[1] = variables[1]
    yCalc = get_sim_edgespectrum_local1(para, postedgearea)
    return (yCalc - y[0])/y[1]


def get_para(inpfile):
    para = []
    for il, line in enumerate(open(inpfile)):
        if "?" in line:
            para.append(float(line.split()[0][1:]))
        else:
            para.append(float(line.split()[0]))
    return np.array(para)


def setFacilityParameter(para):
    sigm1 = para[5] * 10.
    sigm2 = para[6] * 10.
    alph = para[7] * 0.01
    beta = para[8] * 0.01
    sigm = np.sqrt(sigm1**2 + sigm2**2)
    return sigm, alph, beta


#test_fit1234('edge_3.inp.phantom', 2)
#test_fit123('edge_3.inp.phantom', 2)
test_fit('edge_3.inp.phantom')
