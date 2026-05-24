#!/usr/bin/env python
import numpy as np
import textfile
import os
import sys
import pickle
import matplotlib.pyplot as plt
import copy
import datetime
sys.path.append("/home/kazu/ktpro")
from gp_nrca import draw_sample, draw_sample2d, draw_sample2d_mpi
from rits_fit_kt import get_sim_spectrum, get_sim_edgespectrum, get_sim_edgespectrum_local
sys.path.append("/home/kazu/denoise")
import bi2d
np.random.seed(129)


def get_params_edge(rpkl=None, wpkl=None):
    if rpkl:
        with open(rpkl, 'rb') as f:
            param_sets = pickle.load(f)
    else:
        param_sets = {}
        param_sets['param_name'] = ["Tend"]
        param_sets['string'] = ['Te']
        #param_sets['lims_mean'] = [[0.7, 1.0]]    this is set in param_sets_bccrev2_edge.pkl
        param_sets['lims_mean'] = [[0.86062, 0.860621]]
        #param_sets['lims_scale'] = [[0.001, 0.03]] this is set in param_sets_bccrev2_edge.pkl
        param_sets['lims_scale'] = [[0.000001, 0.00003]]
        param_sets['lims_xlim'] = [[10., 30.]]
        param_sets['param_name'].append("Alpha")
        param_sets['string'].append("alpha")
        #param_sets['lims_mean'].append([0.9875, 1.0125]) this is set in param_sets_bccrev2_edge.pkl
        param_sets['lims_mean'].append([1.01521, 1.015211])
        #param_sets['lims_scale'].append([0.001, 0.0025]) this is set in param_sets_bccrev2_edge.pkl
        param_sets['lims_scale'].append([0.000001, 0.0000025])
        param_sets['lims_xlim'].append([10., 30.])
        param_sets['param_name'].append("Beta")
        param_sets['string'].append("beta")
        param_sets['lims_mean'].append([0.4, 1.0])
        param_sets['lims_scale'].append([0.01, 0.04])
        param_sets['lims_xlim'].append([10., 30.])
        param_sets['param_name'].append("Gamma")
        param_sets['string'].append("gamma")
        param_sets['lims_mean'].append([0.9875, 1.0125])
        param_sets['lims_scale'].append([0.001, 0.0025])
        param_sets['lims_xlim'].append([10., 30.])
        param_sets['param_name'].append("dhkl")
        param_sets['string'].append("DHKL")
        #param_sets['lims_mean'].append([2.025, 2.038]) this is set in param_sets_bccrev2_edge.pkl
        param_sets['lims_mean'].append([2.025-0.063, 2.038+0.02])
        param_sets['lims_scale'].append([0.0, 0.000035])
        param_sets['lims_xlim'].append([5., 25.])
        param_sets['param_name'].append("sigma_0")
        param_sets['string'].append("SIGMA_0")
        param_sets['lims_mean'].append([5., 35.])
        param_sets['lims_scale'].append([7.0, 28.0])
        param_sets['lims_xlim'].append([10., 30.])
        param_sets['param_name'].append("mask")
        param_sets['string'].append("MASK")
        param_sets['lims_mean'].append([0.51, 0.52])
        param_sets['lims_scale'].append([0.01, 0.03])
        param_sets['lims_xlim'].append([5., 10.])
        param_sets['maskparams'] = [0.9, 2, 4]
        if wpkl:
            with open(wpkl, 'wb') as f:
                pickle.dump(param_sets, f, 4)
    return param_sets


def draw_params_edge_mpi(rg, param_sets, dim=1, kerneltype='square'):
    if kerneltype == 'random':
        #kerneltypes = ['square', 'single', 'third', 'fifth']
        kerneltypes = ['square', 'third', 'fifth']
        _kerneltype = kerneltypes[rg.integers(3)]
    else:
        _kerneltype = kerneltype
    for pidx in range(len(param_sets['param_name'])):
        lims_mean = param_sets['lims_mean'][pidx]
        lims_scale = param_sets['lims_scale'][pidx]
        lims_xlim = param_sets['lims_xlim'][pidx]
        mean = rg.uniform(low=lims_mean[0], high=lims_mean[1])
        scale = rg.uniform(low=lims_scale[0], high=lims_scale[1])
        xlim = rg.uniform(low=lims_xlim[0], high=lims_xlim[1])
        if dim == 1:
            _params = draw_sample(mean=mean, numsample=1, scale=scale,
                                  xlim=xlim)
        elif dim == 2:
            _params = draw_sample2d_mpi(
                    rg, mean=mean, numsample=1, scale=scale, xlim=xlim,
                    xsize=72, ysize=192, rsize=1, kerneltype=_kerneltype)
        if pidx == 5 and np.min(_params) < 0.: # sigma0
            _params -= np.min(_params)
        if pidx == 0 and np.min(_params) < 0.: # Te
            _params -= (np.min(_params) - 0.05)
        if pidx == 0 and np.max(_params) > 1.:
            _params /= np.max(_params)
        if pidx == 1 and np.min(_params) < 0.: # alpha
            _params -= (np.min(_params) - 0.05)
        #if pidx == 1 and np.max(_params) > 1.:
        #    _params /= np.max(_params)
        if pidx == 2 and np.min(_params) < 0.: # beta
            _params -= (np.min(_params) - 0.05)
        if pidx == 2 and np.max(_params) > 1.:
            _params /= np.max(_params)
        if pidx == 3 and np.min(_params) < 0.: # gamma
            _params -= (np.min(_params) - 0.05)
        #if pidx == 3 and np.max(_params) > 1.:
        #    _params /= np.max(_params)
        if pidx == 6 and np.min(_params) < 0.: # mask
            _params -= np.min(_params)
        if pidx == 6 and np.max(_params) > 1.:
            _params /= np.max(_params)
        if pidx == 0:
            _param_sets = copy.deepcopy(param_sets)
            _param_sets['mean'] = [mean]
            _param_sets['scale'] = [scale]
            _param_sets['xlim'] = [xlim]
            _param_sets['params'] = [_params[0]]
            _param_sets['kerneltype'] = [_kerneltype]
        else:
            _param_sets['mean'].append(mean)
            _param_sets['scale'].append(scale)
            _param_sets['xlim'].append(xlim)
            _param_sets['params'].append(_params[0])
            _param_sets['kerneltype'].append(_kerneltype)
    return _param_sets


def cycles_edge_mpi(pklfile, ns=10, dim=2, kerneltype='square',
                    rpklfile=None,
                    wpklfile=None,
                    ):
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    psize = comm.Get_size()
    import gc
    # for restart job
    if os.path.isfile(pklfile):
        if rank == 0:
            with open(pklfile, 'rb') as f:
                data = pickle.load(f)
                rgs = pickle.load(f)
            data = None
            del data
            gc.collect()
        else:
            rgs = None
        comm.barrier()
        rg = comm.scatter(rgs, root=0)
    else:
        from numpy.random import Generator, PCG64, SeedSequence
        sg = SeedSequence(1234)
        ss = sg.spawn(psize)
        rg = Generator(PCG64(ss[rank]))
    _param_sets_sets = []
    if rpklfile:
        param_sets = get_params_edge(rpkl=rpklfile)
    elif wpklfile:
        param_sets = get_params_edge(wpkl=wpklfile)
    else:
        sys.exit("Error: No rpklfile nor wpklfile are set")
    for ins in range(rank*(ns//psize), (rank+1)*(ns//psize)):
        print('rank=', rank, 'ins=', ins)
        _param_sets_sets.append(draw_params_edge_mpi(rg, param_sets, dim=dim,
                                                     kerneltype=kerneltype))
    comm.barrier()
    param_sets_sets = comm.gather(_param_sets_sets, root=0)
    rg_sets = comm.gather(rg, root=0)
    if rank == 0:
        __param_sets_sets = [__cont for _cont in param_sets_sets for
                             __cont in _cont]
        __param_sets_sets = reform_param_sets(__param_sets_sets)
        if os.path.isfile(pklfile):
            with open(pklfile, 'rb') as f:
                data = pickle.load(f)
                rgs = pickle.load(f)
            __param_sets_sets = data + __param_sets_sets
        if dim == 1:
            with open(pklfile, 'wb') as f:
                pickle.dump(param_sets_sets, f, 4)
        elif dim == 2:
            with open(pklfile, 'wb') as f:
                pickle.dump(__param_sets_sets, f, 4)
                pickle.dump(rg_sets, f, 4)
        print(datetime.datetime.now(), pklfile, ' is saved')


def cycles_edge_mpi_div(nss=10, ns=10, dim=2, kerneltype='square',
                        orgpklfile='param_sets_sets_bccrev2_2d_edge.pkl'):
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    psize = comm.Get_size()
    import gc
    # for restart job
    param_sets = get_params(rpkl='param_sets_bccrev2_edge.pkl')
    for inss in range(nss):
        pklfile = orgpklfile + "." + str(inss+1)
        outpklfile = orgpklfile + "." + str(inss+2)
        #if os.path.isfile(pklfile):
        if rank == 0:
            with open(pklfile, 'rb') as f:
                data = pickle.load(f)
                rgs = pickle.load(f)
            data = None
            del data
            gc.collect()
        else:
            rgs = None
        comm.barrier()
        rg = comm.scatter(rgs, root=0)
        _param_sets_sets = []
        for ins in range(rank*(ns//psize), (rank+1)*(ns//psize)):
            print('rank=', rank, 'ins=', ins)
            _param_sets_sets.append(draw_params_edge_mpi(rg, param_sets,
                                    dim=dim, kerneltype=kerneltype))
        comm.barrier()
        param_sets_sets = comm.gather(_param_sets_sets, root=0)
        rg_sets = comm.gather(rg, root=0)
        if rank == 0:
            __param_sets_sets = [__cont for _cont in param_sets_sets for
                                 __cont in _cont]
            __param_sets_sets = reform_param_sets(__param_sets_sets)
            if dim == 1:
                with open(outpklfile, 'wb') as f:
                    pickle.dump(param_sets_sets, f, 4)
            elif dim == 2:
                with open(outpklfile, 'wb') as f:
                    pickle.dump(__param_sets_sets, f, 4)
                    pickle.dump(rg_sets, f, 4)
            print(datetime.datetime.now(), outpklfile, ' is saved')
        comm.barrier()


def load_param_sets_sets(param_sets_sets_file='param_sets_sets.pkl'):
    with open(param_sets_sets_file, 'rb') as f:
        param_sets_sets = pickle.load(f)
    return param_sets_sets


def run_edge(params, strings, inpfile='edge_3.inp'):
    # Edge version of run_rits.
    for pidx in range(params[0].size):
        os.system('cp edge_3.inp.temp '+inpfile)
        for sidx, string in enumerate(strings):
            _param = "{:.4f}".format(params[sidx].flatten()[pidx])
            textfile.replace(inpfile, string, _param)
        x, y = get_sim_edgespectrum(inpfile=inpfile)
        if pidx == 0:
            bi2d_true = np.zeros((params[0].size, x.shape[0]))
        bi2d_true[pidx] = y
    if params[0].ndim == 2:
        bi2d_true = bi2d_true.reshape((params[0].shape[1], params[0].shape[0],
                                       x.shape[0]))
    return bi2d_true, x


def run_edge_phantom(paramimage, inpfile='edge_3.inp.phantom'):
    # This is for recalculating the rits edge spectra from fitted params.
    # paramimage[paramidx, posidx1, posidx2]]
    counter = 0
    shape = paramimage.shape
    for p1idx in range(shape[1]):
        for p2idx in range(shape[2]):
            if np.sum(np.abs(paramimage[:, p1idx, p2idx])) > 0.:
                counter += 1
                params = paramimage[:, p1idx, p2idx]
                os.system('cp edge_3.inp.temp '+inpfile)
                for sidx, string in enumerate(["A0", "B0", "AHKL", "BHKL",
                                               "DHKL", "SIGMA_0"]):
                    _param = "{:.4f}".format(params[sidx])
                    textfile.replace(inpfile, string, _param)
                x, y = get_sim_edgespectrum(inpfile=inpfile)
                if counter == 1:
                    bi2d_true = np.ones((shape[1], shape[2], x.shape[0]))
                if y.shape[0] == bi2d_true.shape[2]:
                    bi2d_true[p1idx, p2idx] = y
    return bi2d_true, x


def run_edge_local_phantom(paramimage, inpfile='edge_3.inp.phantom'):
    # This is for recalculating the rits edge spectra from fitted params.
    # paramimage[paramidx, posidx1, posidx2]]
    # Simulated spectra are obtained by python methods wo the rits code.
    counter = 0
    shape = paramimage.shape
    for p1idx in range(shape[1]):
        for p2idx in range(shape[2]):
            if np.sum(np.abs(paramimage[:, p1idx, p2idx])) > 0.:
                counter += 1
                params = paramimage[:, p1idx, p2idx]
                os.system('cp edge_3.inp.temp '+inpfile)
                for sidx, string in enumerate(["A0", "B0", "AHKL", "BHKL",
                                               "DHKL", "SIGMA_0"]):
                    _param = "{:.4f}".format(params[sidx])
                    textfile.replace(inpfile, string, _param)
                x, y = get_sim_edgespectrum_local(inpfile=inpfile)
                if counter == 1:
                    bi2d_true = np.ones((shape[1], shape[2], y.shape[0]))
                if y.shape[0] == bi2d_true.shape[2]:
                    bi2d_true[p1idx, p2idx] = y
    return bi2d_true


def run_edge_local2_phantom(paramimage, inpfile='edge_3.inp.phantom'):
    # This is for recalculating the rits edge spectra from fitted params.
    # paramimage[paramidx, posidx1, posidx2]]
    # Simulated spectra are obtained by python methods wo the rits code.
    para = get_para(inpfile)
    counter = 0
    shape = paramimage.shape
    for p1idx in range(shape[1]):
        for p2idx in range(shape[2]):
            if np.sum(np.abs(paramimage[:, p1idx, p2idx])) > 0.:
                counter += 1
                para[0:6] = paramimage[0:6, p1idx, p2idx]
                y = get_sim_edgespectrum_local(para)
                if counter == 1:
                    bi2d_true = np.ones((shape[1], shape[2], y.shape[0]))
                #if y.shape[0] == bi2d_true.shape[2]:
                #    bi2d_true[p1idx, p2idx] = y
                if np.sum(np.isnan(y)) > 0:
                    bi2d_true[p1idx, p2idx] = np.ones_like(y) 
                else:
                    bi2d_true[p1idx, p2idx] = y
    return bi2d_true


def run_edge_local2_phantom_mean(paramimage, inpfile='edge_3.inp.phantom'):
    para = get_para(inpfile)
    counter = 0
    for i in range(4):
        paramimage[i][paramimage[i]!=0] = np.round(np.mean(paramimage[i][paramimage[i]!=0]), 2)
    shape = paramimage.shape
    for p1idx in range(shape[1]):
        for p2idx in range(shape[2]):
            if np.sum(np.abs(paramimage[:, p1idx, p2idx])) > 0.:
                counter += 1
                para[0:6] = paramimage[0:6, p1idx, p2idx]
                x, y = get_sim_edgespectrum_local(para)
                if counter == 1:
                    bi2d_true = np.ones((shape[1], shape[2], x.shape[0]))
                if y.shape[0] == bi2d_true.shape[2]:
                    bi2d_true[p1idx, p2idx] = y
    return bi2d_true, x


def get_para(inpfile):
    para = []
    for il, line in enumerate(open(inpfile)):
        if "?" in line:
            para.append(float(line.split()[0][1:]))
        else:
            para.append(float(line.split()[0]))
    return np.array(para)


def single_edgecomputation(pklfile='param_sets_sets_bccrev_edge.pkl'):
    param_sets_sets = load_param_sets_sets(param_sets_sets_file=pklfile)
    #param_sets_sets = reform_param_sets(param_sets_sets)
    iniidx = int(sys.argv[1])
    param_sets = param_sets_sets[iniidx]
    inpfile = 'edge_3.inp.' + str(iniidx)
    # param_sets are used except for mask
    bi3d_true, x = run_edge(param_sets['params'][:-1],  param_sets['string'][:-1],
                            inpfile=inpfile)
    # apply mask
    # (192,72,152) (72,192)
    mask = param_sets['params'][-1]
    shape = bi3d_true.shape
    bi3d_true = bi3d_true.reshape(shape[1], shape[0], shape[2]).transpose((2, 0, 1))
    bi3d_true = np.ones_like(bi3d_true)*(1. - mask) + bi3d_true * mask
    bi3d_true = bi3d_true.transpose((1, 2, 0))
    with open('bi2dedge.pkl.'+str(iniidx), 'wb') as f:
        pickle.dump(bi3d_true, f, 4)
        pickle.dump(x, f, 4)


def reform_param_sets(param_sets_sets):
    from scipy.ndimage import gaussian_filter
    _params = []
    for param_sets in param_sets_sets:
        _params.append(param_sets['params'])
    _params = np.array(_params)
    Te = _params[:, 0]
    alpha = _params[:, 1]
    beta = _params[:, 2]
    gamma = _params[:, 3]
    mask = _params[:, 6]
    maskparams = param_sets_sets[0]['maskparams']
    print("CHK", mask.shape)
    for didx in range(mask.shape[0]):
        ave = np.average(mask[didx])
        mask[didx][mask[didx] > ave*maskparams[0]] = 1.
        mask[didx][mask[didx] <= ave*maskparams[0]] = 0.
        mask[didx] = gaussian_filter(
                mask[didx], sigma=np.random.uniform(low=maskparams[1],
                                                    high=maskparams[2]))

    Tiext = alpha*Te
    Ti = beta*Tiext
    Teext = gamma*Ti
    tini = 23000.
    tend = 26020.
    print('chk', np.min(Ti), np.min(Te), np.min(Tiext), np.min(Teext))
    a0 = (np.log(Te/Tiext)*tini + (tini - tend)*np.log(Tiext)) / (tend - tini)
    b0 = 10000.*np.log(Te/Tiext) / (tini - tend)
    ahkl = (tini*np.log(Teext*Tiext/Ti/Te) + (tend - tini)*np.log(Tiext/Ti)) / (tend - tini)
    bhkl = 10000.*np.log(Teext*Tiext/Ti/Te) / (tini - tend)
    for idx in range(a0.shape[0]):
        param_sets_sets[idx]['params'][0] = a0[idx]
        param_sets_sets[idx]['params'][1] = b0[idx]
        param_sets_sets[idx]['params'][2] = ahkl[idx]
        param_sets_sets[idx]['params'][3] = bhkl[idx]
        param_sets_sets[idx]['params'][6] = mask[idx]
        param_sets_sets[idx]['string'][0] = "A0"
        param_sets_sets[idx]['string'][1] = "B0"
        param_sets_sets[idx]['string'][2] = "AHKL"
        param_sets_sets[idx]['string'][3] = "BHKL"
        param_sets_sets[idx]['param_name'][0] = "a0"
        param_sets_sets[idx]['param_name'][1] = "b0"
        param_sets_sets[idx]['param_name'][2] = "ahkl"
        param_sets_sets[idx]['param_name'][3] = "bhkl"
    return param_sets_sets


def single_edgecomputation_phantom(pklfile='paramimage_expt_5models.pkl'):
    # for constructing phantom data from the result of the rits edge fitting.
    with open(pklfile, 'rb') as f:
        paramimage = pickle.load(f)
    inpfile = 'edge_3.inp.phantom'
    #bi3d_true, x = run_edge_local2_phantom_mean(paramimage, inpfile=inpfile)
    #import matplotlib.pyplot as plt
    #plt.imshow(np.sum(np.abs(paramimage), axis=0))
    #plt.show()
    bi3d_true = run_edge_local2_phantom(paramimage, inpfile=inpfile)
    #bi3d_true, x = run_edge_phantom(paramimage, inpfile=inpfile)
    #with open('bi2dsingle_scratch_rev.pkl.phantom_local', 'wb') as f:
    #with open('bi2dsingle_denoised_5models.pkl.phantom_local', 'wb') as f:
    #with open('bi2dsingle_scratch_rev2.pkl.phantom_local', 'wb') as f:
    #with open('bi2dsingle_denoised_rev2_5models.pkl.phantom_local', 'wb') as f:
    #with open('bi2dsingle_scratch_rev3.pkl.phantom_local', 'wb') as f:
    with open('bi2dsingle_denoised_5models_ktrand.pkl.phantom_local', 'wb') as f:
        #with open('bi2dsingle_scratch_rev4.pkl.phantom_local', 'wb') as f:
        #with open('bi2dsingle_scratch.pkl.phantom', 'wb') as f:
        pickle.dump(bi3d_true, f, 4)


def synthesize_bi3ddata():
    with open('/home/kazu/desktop/240424/uNID_data_KO/211/openbeam.pkl',
              'rb') as f:
        openbeamexp = pickle.load(f)
        openbeamexp_noisy = pickle.load(f)
        sampleexp = pickle.load(f)
        sampleexp_noisy = pickle.load(f)
    with open('bi3d_test_single_rev2_true_edge_ktrand_wosingle.pkl', 'rb') as f:
        datasets = pickle.load(f)
    #data = datasets['target'].transpose((0,2,1,3))
    #shape = datasets['target'].shape
    #data = datasets['target'].reshape((shape[0], shape[2], shape[1], shape[3]))
    data = datasets['target']
    #plt.imshow(data[0,:,:,:].sum(axis=2))
    #plt.show()
    print("openbeamexp:", openbeamexp.shape)
    print("data:", data.shape)
    _shape = openbeamexp.shape
    # open beam intensities integrated w.r.t tof
    tcount_openbeamexp = openbeamexp.sum(axis=0)
    # noisy open beam intensities integrated w.r.t tof
    tcount_openbeamexp_noisy = openbeamexp_noisy.sum(axis=0)
    # spatially averaged open beam spectrum
    mean_openbeamexp = np.mean(openbeamexp, axis=(1, 2))*1.15
    sample = (data*mean_openbeamexp).transpose((0, 3, 2, 1))\
        * tcount_openbeamexp / np.sum(mean_openbeamexp)*1.5
    op = ((np.zeros_like(data[0])+1)*mean_openbeamexp).transpose((2, 1, 0))\
        * tcount_openbeamexp/np.sum(mean_openbeamexp)
    print("sample:", sample.shape)
    print("op:", op.shape)
    x = np.arange(_shape[0])*20.+2300
    output = np.vstack((x, sample[10, :, 40, 12]/op[:, 40, 12]))
    output = np.vstack((output, np.zeros_like(x) + 0.01))
    np.savetxt('check_simuspectrum.txt', output.T)

    px = 57
    py = 44
    snum = 4
    # for median filter
    k = 24
    w = sampleexp.shape[0]
    idx = np.fromfunction(lambda i, j: i + j, (k, w), dtype=np.int64) - k // 2
    idx[idx < 0] = 0
    idx[idx > w - 1] = w - 1
    #
    smoothed = np.median(sampleexp[idx, 40-8, 24-16], axis=0)
    fig, ax = plt.subplots(2, 2)
    for didx, (spectrum, label) in enumerate(zip([sampleexp[:, px, py],
                                             np.random.poisson(
                                                 sample[snum, :, px, py])],
                                             ['sampleexp', 'samplesim'])):
        smoothed = np.median(spectrum[idx], axis=0)
        ax[0, didx].plot(spectrum, label=label)
        ax[0, didx].plot(smoothed, label='medianfiltered')
        ax[0, didx].set_xlabel(r'$\lambda$ / ch')
        ax[0, didx].set_xlabel('TOF / ch')
        ax[0, didx].set_ylabel('Neutron Count')
        ax[0, didx].legend()
        ax[1, didx].plot(spectrum - smoothed)
        ax[1, didx].set_ylabel(label + '-medianfiltered')
        ax[1, didx].set_xlabel('TOF / ch')
    plt.tight_layout()
    #plt.savefig('noise_levels_sample_exp_simu.png')
    plt.show()

    fig, ax = plt.subplots(4, 2, figsize=(10, 15))
    for didx, (spectrum, label) in enumerate(zip([sampleexp[:, px, py],
                                                  sampleexp_noisy[:, px, py],
                                                  openbeamexp[:, px, py],
                                                  openbeamexp_noisy[:, px, py]
                                                  ],
                                                 ['sample',
                                                  'sample_part',
                                                  'openbeam',
                                                  'openbeam_part']
                                                 )):
        smoothed = np.median(spectrum[idx], axis=0)
        ax[didx, 0].plot(spectrum, label=label)
        ax[didx, 0].plot(smoothed, label='medianfiltered')
        ax[didx, 0].set_xlabel(r'$\lambda$ / ch')
        ax[didx, 0].set_xlabel('TOF / ch')
        ax[didx, 0].set_ylabel('Neutron Count')
        ax[didx, 0].legend()
        ax[didx, 1].plot(spectrum - smoothed)
        ax[didx, 1].set_ylabel('Diff.')
    plt.tight_layout()
    #plt.savefig('noise_levels_sample_exp_full_part.png')
    plt.show()

    fig, ax = plt.subplots(4, 2, figsize=(10, 10))
    _op = np.random.poisson(op)
    vmin1 = np.min(np.vstack((_op[:, px, :], openbeamexp[:, px, :])))
    vmax1 = np.max(np.vstack((_op[:, px, :], openbeamexp[:, px, :])))
    vmin2 = np.min(np.vstack((_op[:, :, py], openbeamexp[:, :, py])))
    vmax2 = np.max(np.vstack((_op[:, :, py], openbeamexp[:, :, py])))
    vmin3 = np.min(np.vstack((_op.sum(axis=0), openbeamexp.sum(axis=0))))
    vmax3 = np.max(np.vstack((_op.sum(axis=0), openbeamexp.sum(axis=0))))
    ymax = np.max(np.vstack((_op[:, px, py], openbeamexp[:, px, py])))
    ymin = np.min(np.vstack((_op[:, px, py], openbeamexp[:, px, py])))
    for didx, (image, name) in enumerate(zip([_op, openbeamexp],
                                             ['simulated_openbeam',
                                                 'expt_openbeam'])):
        ax[0, didx].set_title(name + f' crosssection at y={px} ch')
        ax[0, didx].set_xlabel(r'$\lambda$ / ch')
        ax[0, didx].set_ylabel('y / ch')
        ax[0, didx].imshow(image[:, px, :].T, vmin=vmin1, vmax=vmax1, aspect=5)
        ax[1, didx].set_title(name + f' crosssection at x={py} ch')
        ax[1, didx].set_xlabel(r'$\lambda$ / ch')
        ax[1, didx].imshow(image[:, :, py].T, vmin=vmin2, vmax=vmax2, aspect=10)
        ax[1, didx].set_ylabel('x / ch')
        ax[2, didx].set_title(name + ' projected onto xy')
        ax[2, didx].imshow(image.sum(axis=0), vmin=vmin3, vmax=vmax3)
        ax[2, didx].set_xlabel('x / ch')
        ax[2, didx].set_ylabel('y / ch')
        ax[2, didx].axhline(px, color='r', linestyle='--', lw=1)
        ax[2, didx].axvline(py, color='r', linestyle='--', lw=1)
        ax[3, didx].set_title(name + f' intensity at y={px} ch, x={py} ch')
        ax[3, didx].plot(image[:, px, py])
        ax[3, didx].set_ylim([ymin, ymax])
        ax[3, didx].set_xlabel(r'$\lambda$ / ch')
        ax[3, didx].set_ylabel('Neutron Count')
    plt.tight_layout()
    #plt.savefig('openbeam_data_expt_simu.png')
    plt.show()
    fig, ax = plt.subplots(4, 2, figsize=(10, 10))
    _sample = np.random.poisson(sample[snum])
    vmin1 = np.min(np.vstack((_sample[:, px, :], sampleexp[:, px, :])))
    vmax1 = np.max(np.vstack((_sample[:, px, :], sampleexp[:, px, :])))
    vmin2 = np.min(np.vstack((_sample[:, :, py], sampleexp[:, :, py])))
    vmax2 = np.max(np.vstack((_sample[:, :, py], sampleexp[:, :, py])))
    vmin3 = np.min(np.vstack((_sample.sum(axis=0), sampleexp.sum(axis=0))))
    vmax3 = np.max(np.vstack((_sample.sum(axis=0), sampleexp.sum(axis=0))))
    ymax = np.max(np.vstack((_sample[:, px, py], sampleexp[:, px, py])))
    ymin = np.min(np.vstack((_sample[:, px, py], sampleexp[:, px, py])))
    for didx, (image, name) in enumerate(zip([_sample, sampleexp],
                                             ['simulated_sample',
                                                 'expt_sample'])):
        ax[0, didx].set_title(name + ' crosssection at y=4 ch')
        ax[0, didx].set_xlabel(r'$\lambda$ / ch')
        ax[0, didx].set_ylabel('y / ch')
        ax[0, didx].imshow(image[:, px, :].T, vmin=vmin1, vmax=vmax1, aspect=5)
        ax[1, didx].set_title(name + ' crosssection at x=32 ch')
        ax[1, didx].set_xlabel(r'$\lambda$ / ch')
        ax[1, didx].imshow(image[:, :, py].T, vmin=vmin2, vmax=vmax2, aspect=10)
        ax[1, didx].set_ylabel('x / ch')
        ax[2, didx].set_title(name + ' projected onto xy')
        ax[2, didx].imshow(image.sum(axis=0), vmin=vmin3, vmax=vmax3)
        ax[2, didx].set_xlabel('x / ch')
        ax[2, didx].set_ylabel('y / ch')
        ax[2, didx].axhline(px, color='r', linestyle='--', lw=1)
        ax[2, didx].axvline(py, color='r', linestyle='--', lw=1)
        ax[3, didx].set_title(name + ' intensity at y=4 ch, x=32 ch')
        ax[3, didx].plot(image[:, px, py])
        ax[3, didx].set_ylim([ymin, ymax])
        ax[3, didx].set_xlabel(r'$\lambda$ / ch')
        ax[3, didx].set_ylabel('Neutron Count')
    plt.tight_layout()
    plt.savefig('sample_data_expt_simu.png')
    plt.show()
    sample = sample.transpose((0, 3, 2, 1))[:, np.newaxis].astype('float32')
    #op = op.transpose((2, 1, 0))
    print(sample.shape)
    print(data.shape)
    #_shape2 = sample.shape
    #sample[:, 1::2] = np.flip(sample[:, 1::2], axis=2)
    #data[:, 1::2] = np.flip(data[:, 1::2], axis=2)
    #sample = sample.reshape((-1, _shape2[2]*2, _shape2[3]))
    #data = data.reshape((-1, _shape2[2]*2, _shape2[3]))
    #fig, ax = plt.subplots(4, 1, figsize=(10, 10))
    #ax[0].imshow(sample[0])
    #ax[1].imshow(data[0])
    #ax[2].plot(data[0, :, 300], marker='o')
    #ax[3].plot(sample[0, :, 300], marker='x')
    #plt.show()
    #sim_datasets = {}
    #sim_datasets['sample'] = sample
    #sim_datasets['op'] = op
    #sim_datasets['x'] = x
    #with open('bi3d_testbcc_simudata_rev2_lim_single_resize_full_211_true_edgev_ktrand_wosingle.pkl.12', 'wb') as f:
    #    pickle.dump(sample, f, 4)


def synthesize_bi3ddata_with_various_intensities():
    with open('/home/kazu/desktop/240424/uNID_data_KO/211/openbeam.pkl',
              'rb') as f:
        openbeamexp = pickle.load(f)
    with open('/home/kazu/desktop/240424/connect2d/sampled/bi3d_test_single_rev2_true_edge_ktrand.pkl', 'rb') as f:
        data = pickle.load(f)['target']
    print(data.shape)
    # open beam intensities integrated w.r.t tof
    tcount_openbeamexp = openbeamexp.sum(axis=0)
    # noisy open beam intensities integrated w.r.t tof
    # spatially averaged open beam spectrum
    #for amp in 2**np.linspace(1, 10, 10):
    for amp in [2., 4., 8., 16.]:
        mean_openbeamexp = np.mean(openbeamexp, axis=(1, 2))*1.15
        sample = (data*mean_openbeamexp).transpose((0, 3, 2, 1))\
            * tcount_openbeamexp / np.sum(mean_openbeamexp)/315715.*553690.*amp
        sample = sample.transpose((0, 3, 2, 1))[:, np.newaxis].astype('float32')
        op = ((np.zeros_like(data[0])+1)*mean_openbeamexp).transpose((2, 1, 0))\
            * tcount_openbeamexp/np.sum(mean_openbeamexp) * amp
        with open('bi3d_testbcc_simudata_rev2_lim_single_resize_full_211_true_edgev_ktrand_intx' + str(amp)[:-2] + '.pkl', 'wb') as f:
            pickle.dump(sample, f, 4)


def show_histograms_of_rits_params(param_sets_sets_file, paramname):
    with open(param_sets_sets_file, 'rb') as f:
        params = pickle.load(f)
    print("number of param sets:", len(params))
    paramvals = []
    pidx = 0
    print(params[0].keys())
    print(params[0]['lims_mean'])
    #plt.imshow(params[0]['params'][4])
    #plt.show()
    for pidx, name in enumerate(params[0]['param_name']):
        if name == paramname:
            target_pidx = pidx
    print('CHECK', target_pidx)
    for pset in params:
        paramvals.append(pset['params'][target_pidx])
    paramvals = np.asarray(paramvals)
    plt.hist(paramvals.flatten(), bins=40)
    plt.show()


def calculate_alpha_and_Te_for_optimum_a0_b0(a0=0.025, b0=0.075):
    tini = 23000.
    tend = 26020.
    Te = np.exp(-a0-0.0001*b0*tend)
    Tiext = np.exp(-a0-0.0001*b0*tini)
    alpha = Tiext/Te
    print('alpha=', alpha)
    print('Te=', Te)




# cycles_edge_mpi(ns=1280, dim=2, kerneltype='random', pklfile='param_sets_sets_bccrev3_2d_single_edge_true_edge_ktrand_wosingle.pkl')
# cycles_edge_mpi_div(nss=11, ns=1280, dim=2, kerneltype='random', orgpklfile='param_sets_sets_bccrev2_2d_single_edge_true_edge_ktrand_wosingle.pkl')
# single_edgecomputation(pklfile='param_sets_sets_bccrev2_2d_single_edge_true_edge.pkl')
# single_edgecomputation_phantom()
# synthesize_bi3ddata()
# check_data()
# crude_parallel_computation()
