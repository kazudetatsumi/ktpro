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
#from gp_nrca import draw_sample, draw_sample2d, draw_sample2d_mpi
from gp_nrca import draw_sample, draw_sample2d_mpi
from rits_fit_kt import get_sim_spectrum
sys.path.append("/home/kazu/denoise")
import bi2d
np.random.seed(129)


def get_params():
    param_sets = {}
    param_sets['param_name'] = ["lattice constance"]
    param_sets['string'] = ['LATCON']
    ##param_sets['lims_mean'] = [[2.5, 3.2]]
    #param_sets['lims_mean'] = [[2.86338, 2.86340]]
    #param_sets['lims_mean'] = [[2.85315, 2.908040]]
    #param_sets['lims_mean'] = [[2.86, 2.87]]
    param_sets['lims_mean'] = [[2.86, 2.88]]
    ##param_sets['lims_scale'] = [[0.0, 0.04]]
    #param_sets['lims_scale'] = [[0.0, 0.00005]]
    param_sets['lims_scale'] = [[0.0, 0.00006]]
    ##param_sets['lims_scale'] = [[0.0, 0.00001]]
    #param_sets['lims_xlim'] = [[10., 30.]]
    param_sets['lims_xlim'] = [[5., 25.]]
    #45 Crsytallite size, 2.0〜9.0 ミクロン. エッジのジャンプ高さが低くなる。
    param_sets['param_name'].append("crsytallite size")
    param_sets['string'].append("CRSIZE")
    ##param_sets['lims_mean'].append([0.5, 1.1])
    #param_sets['lims_mean'].append([0.800057, 0.800059])
    #param_sets['lims_mean'].append([0.186038, 1.31959])
    param_sets['lims_mean'].append([0.5, 2.5])
    ##param_sets['lims_scale'].append([0.1, 2.0])
    ##param_sets['lims_scale'].append([0.0, 0.00001])
    #param_sets['lims_scale'].append([0.0, 0.01])
    param_sets['lims_scale'].append([0.0, 0.10])
    #param_sets['lims_xlim'].append([10., 30.])
    param_sets['lims_xlim'].append([20., 40.])
    #23 Microstrain in Jorgensen function, エッジの傾きが変わる。
    # 0〜6ぐらい．
    param_sets['param_name'].append("microstrain in Jorgensen func")
    param_sets['string'].append("MICRST")
    ##param_sets['lims_mean'].append([0.0, 6.0])
    #param_sets['lims_mean'].append([1.65100, 1.65102])
    #param_sets['lims_mean'].append([1.0, 15.7124])
    param_sets['lims_mean'].append([2.0, 15.])
    ##param_sets['lims_scale'].append([0.1, 2.0])
    ##param_sets['lims_scale'].append([0.0, 0.00001])
    #param_sets['lims_scale'].append([0.1, 2.0])
    param_sets['lims_scale'].append([0.1, 3.0])
    #param_sets['lims_xlim'].append([10., 30.])
    param_sets['lims_xlim'].append([10., 40.])
    #24 Crysallite size in Jorgensen function, エッジの傾きが変わる。
    # 0〜6ぐらい．
    #param_sets['param_name'].append("crsytallite size in Jorgensen func")
    #param_sets['string'].append("CRSIZJ")
    #param_sets['lims_mean'].append([0.0, 6.0])
    #param_sets['lims_scale'].append([0.1, 2.0])
    #param_sets['lims_xlim'].append([10., 30.])
    #32 March-Dollase 係数,  1.5〜2.0 (蘇さんの歯車), 0.3〜0.9 (佐藤さんの溶接材)。
    # <1でビームに平行に成長、>1で垂直に成長。
    # 特定のhkl列のエッジジャンプが変わる．
    param_sets['param_name'].append("March-Dollase Coefficient")
    param_sets['string'].append("MDCOEF")
    #param_sets['lims_mean'].append([0.3, 2.0])
    ##param_sets['lims_mean'].append([1.5, 2.0])
    #param_sets['lims_mean'].append([1.942830, 1.94283])
    #param_sets['lims_mean'].append([1.33567, 2.53587])
    param_sets['lims_mean'].append([1.5, 4.])
    ##param_sets['lims_scale'].append([0.01, 0.1])
    #param_sets['lims_scale'].append([0.0, 0.000001])
    #param_sets['lims_scale'].append([0.01, 2.0])
    param_sets['lims_scale'].append([0.01, 1.0])
    #param_sets['lims_xlim'].append([10., 30.])
    param_sets['lims_xlim'].append([10., 40.])
    # 元素特性1
    #param_sets['param_name'].append("Coherent Scattering Length")
    #param_sets['string'].append("COHESL")
    ##param_sets['lims_mean'].append([9.3, 9.6])
    #param_sets['lims_mean'].append([9.44, 9.46])
    ##param_sets['lims_scale'].append([0.01, 0.1])
    #param_sets['lims_scale'].append([0.0, 0.000001])
    #param_sets['lims_xlim'].append([10., 30.])
    ## 投影原子数密度
    param_sets['param_name'].append("Projected Atomic Number Density")
    param_sets['string'].append("PRODEN")
    ##param_sets['lims_mean'].append([1.0, 7.0])
    #param_sets['lims_mean'].append([3.50487, 3.50489])
    #param_sets['lims_mean'].append([1.33478, 3.96254])
    param_sets['lims_mean'].append([1., 4.])
    ##param_sets['lims_scale'].append([0.05, 0.5])
    #param_sets['lims_scale'].append([0.0, 0.000005])
    #param_sets['lims_scale'].append([0.05, 2.0])
    param_sets['lims_scale'].append([0.05, 4.0])
    #param_sets['lims_xlim'].append([6., 10.])
    param_sets['lims_xlim'].append([16., 20.])

    return param_sets


def draw_params(param_sets, dim=1):
    for pidx in range(len(param_sets['param_name'])):
        lims_mean = param_sets['lims_mean'][pidx]
        lims_scale = param_sets['lims_scale'][pidx]
        lims_xlim = param_sets['lims_xlim'][pidx]
        mean = np.random.uniform(low=lims_mean[0], high=lims_mean[1])
        scale = np.random.uniform(low=lims_scale[0], high=lims_scale[1])
        xlim = np.random.uniform(low=lims_xlim[0], high=lims_xlim[1])
        if dim == 1:
            _params = draw_sample(mean=mean, numsample=1, scale=scale,
                                  xlim=xlim)
        elif dim == 2:
            _params = draw_sample2d(mean=mean, numsample=1, scale=scale,
                                    xlim=xlim, xsize=96, rsize=96)
        if np.min(_params) < 0.:
            _params -= np.min(_params)
        if pidx == 0:
            _param_sets = copy.deepcopy(param_sets)
            _param_sets['mean'] = [mean]
            _param_sets['scale'] = [scale]
            _param_sets['params'] = [_params[0]]
        else:
            _param_sets['mean'].append(mean)
            _param_sets['scale'].append(scale)
            _param_sets['params'].append(_params[0])
    return _param_sets


def draw_params_mpi(rg, param_sets, dim=1):
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
            _params = draw_sample2d_mpi(rg, mean=mean, numsample=1,
                                    scale=scale, xlim=xlim, xsize=96,
                                    rsize=96)
        #if pidx==0:
        #    print(param_sets['param_name'][pidx], param_sets['lims_mean'][pidx], mean, _params)
        if np.min(_params) < 0.:
            _params -= np.min(_params)
        if pidx == 0:
            _param_sets = copy.deepcopy(param_sets)
            _param_sets['mean'] = [mean]
            _param_sets['scale'] = [scale]
            _param_sets['params'] = [_params[0]]
        else:
            _param_sets['mean'].append(mean)
            _param_sets['scale'].append(scale)
            _param_sets['params'].append(_params[0])
    return _param_sets


def cycles(ns=2, dim=1):
    param_sets_sets = []
    param_sets = get_params()
    for ins in range(ns):
        param_sets_sets.append(draw_params(param_sets, dim=dim))
    if dim == 1:
        with open('param_sets_sets_bccrev3.pkl', 'wb') as f:
            pickle.dump(param_sets_sets, f, 4)
    elif dim == 2:
        with open('param_sets_sets_2d.pkl', 'wb') as f:
            pickle.dump(param_sets_sets, f, 4)


def cycles_mpi(ns=10, dim=2, pklfile='param_sets_sets_2dlarge.pkl'):
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
    param_sets = get_params()
    for ins in range(rank*(ns//psize), (rank+1)*(ns//psize)):
        _param_sets_sets.append(draw_params_mpi(rg, param_sets, dim=dim))
    comm.barrier()
    param_sets_sets = comm.gather(_param_sets_sets, root=0)
    rg_sets = comm.gather(rg, root=0)
    if rank == 0:
        __param_sets_sets = [__cont for _cont in param_sets_sets for
                             __cont in _cont]
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


def divide_paramdata(pklfile='param_sets_sets_2dlarge.pkl'):
    with open(pklfile, 'rb') as f:
        param_sets_sets = pickle.load(f)
        rgs = pickle.load(f)
    for pidx, param_sets in enumerate(param_sets_sets):
        with open("divided_params/"+pklfile+"."+str(pidx), 'wb') as f:
            pickle.dump(param_sets, f, 4)


def load_param_sets_sets(param_sets_sets_file='param_sets_sets.pkl'):
    with open(param_sets_sets_file, 'rb') as f:
        param_sets_sets = pickle.load(f)
    return param_sets_sets


def run_rits(params, strings, inpfile='rits_initial.inp'):
    #for pidx in range(len(params[0])):
    for pidx in range(params[0].size):
        os.system('cp rits_initial.inp.temp '+inpfile)
        for sidx, string in enumerate(strings):
            _param = "{:.4f}".format(params[sidx].flatten()[pidx])
            textfile.replace(inpfile, string, _param)
            if string == "COHESL" and pidx == 0:
                params2 = params[sidx]**2*4*np.pi/100.
            if string == "COHESL":
                _param2 = "{:.4f}".format(params2.flatten()[pidx])
                textfile.replace(inpfile, "COHXSC", _param2)
        x, y = get_sim_spectrum(inpfile=inpfile)
        if pidx == 0:
            bi2d_true = np.zeros((params[0].size, x.shape[0]))
        bi2d_true[pidx] = y
    if params[0].ndim == 2:
        bi2d_true = bi2d_true.reshape((params[0].shape[0], params[0].shape[1], x.shape[0]))
    return bi2d_true, x


def get_noisydata(bi2d, x, timescale):
    return np.random.poisson(bi2d*timescale*x)/x


#cycles(ns=10000)

# routine for crude parallel computations of rits
def crude_parallel_computation(pklfile='param_sets_sets_bccrev2.pkl'):
    inino = int(sys.argv[1])
    inpfile = 'rits_initial.inp.' + str(inino)
    param_sets_sets = load_param_sets_sets(param_sets_sets_file=pklfile)
    for pid, param_sets in enumerate(param_sets_sets[inino*100:inino*100+100]):
        bi2d_true, x = run_rits(param_sets['params'],  param_sets['string'],
                                inpfile=inpfile)
        if pid == 0:
            bi2dt = np.zeros((100, bi2d_true.shape[0], bi2d_true.shape[1]))
        bi2dt[pid] = bi2d_true
    with open('/home/kazu/desktop/240424/bi2d/' + str(inino) + '_rev2_chk.pkl', 'wb') as f:
        pickle.dump(bi2dt, f, 4)
        pickle.dump(x, f, 4)
    #bi2d_noisy = get_noisydata(bi2d_true, x, 200)
    #bi2d(bi2d_true, bi2d_noisy)


# routine for mpi parallel computations of rits
def mpi_parallel_computation(pklfile='param_sets_sets_2dlarge.pkl'):
    print('loading ', pklfile)
    param_sets_sets = load_param_sets_sets(param_sets_sets_file=pklfile)
    print('loaded ', pklfile)
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    psize = comm.Get_size()
    inpfile = 'rits_initial.inp.' + str(rank)
    #ns = len(param_sets_sets)
    iniidx = int(sys.argv[1])
    ns = 1
    nskip = ns*iniidx
    if rank == 0:
        print("NSKIP:", nskip)
    for splidx in range(rank*(ns//psize), (rank+1)*(ns//psize)):
        param_sets = param_sets_sets[splidx+nskip]
        bi2d_true, x = run_rits(param_sets['params'],  param_sets['string'],
                                inpfile=inpfile)
        if splidx == rank*(ns//psize):
            spectrar = np.zeros((ns//psize, bi2d_true.flatten().shape[0]))
        spectrar[splidx - rank*ns//psize, :] = bi2d_true.flatten()
    comm.barrier()
    bi2dt = np.array(comm.gather(spectrar.flatten(), root=0))
    if rank == 0:
        if bi2d_true.ndim == 2:
            bi2dt = bi2dt.reshape((ns, bi2d_true.shape[0], bi2d_true.shape[1]))
        elif bi2d_true.ndim == 3:
            bi2dt = bi2dt.reshape((ns, bi2d_true.shape[0], bi2d_true.shape[1], bi2d_true.shape[2]))
        with open('/home/kazu/desktop/240424/bi2d/bi2dlarge_test_.pkl.'+str(iniidx),
                  'wb') as f:
            pickle.dump(bi2dt, f, 4)
            pickle.dump(x, f, 4)


def single_computation(pklfile='param_sets_sets_bccrev.pkl'):
    iniidx = int(sys.argv[1])
    #pklfile = "divided_params/"+pklfile+"."+str(iniidx)
    #param_sets = load_param_sets_sets(param_sets_sets_file=pklfile)
    param_sets = load_param_sets_sets(param_sets_sets_file=pklfile)[iniidx]
    inpfile = 'rits_initial.inp.' + str(iniidx)
    bi3d_true, x = run_rits(param_sets['params'],  param_sets['string'],
                            inpfile=inpfile)
    with open('/home/kazu/desktop/240424/bi2dsingle.pkl.'+str(iniidx),
              'wb') as f:
        pickle.dump(bi3d_true, f, 4)
        pickle.dump(x, f, 4)


def synthesize_bi2ddata():
    with open('/home/kazu/desktop/240424/uNID_data_KO/433/openbeam.pkl',
              'rb') as f:
        openbeamexp = pickle.load(f)[:, 8:-5, 16:-16]
        openbeamexp_noisy = pickle.load(f)[:, 8:-5, 16:-16]
        sampleexp = pickle.load(f)[:, 8:-5, 16:-16]
        sampleexp_noisy = pickle.load(f)[:, 8:-5, 16:-16]
        print("CHK", openbeamexp.shape)
    with open('/home/kazu/desktop/240424/bi2d/bi2d_testbcc_rev2.pkl',
              'rb') as f:
        datasets = pickle.load(f)
    data = datasets['target'][:, 8:-5, :]
    print(data.shape)
    _shape = openbeamexp.shape
    tcount_openbeamexp = np.tile(openbeamexp.sum(axis=0),
                                 (1, 1612)).T[:data.shape[0]]
    tcount_openbeamexp_noisy = np.tile(openbeamexp_noisy.sum(axis=0),
                                       (1, 1612)).T[:data.shape[0]]
    mean_openbeamexp = openbeamexp.transpose((2, 1, 0))\
        .sum(axis=0).sum(axis=0)/_shape[1]/_shape[2]*1.15
    sample = (data*mean_openbeamexp).transpose((2, 0, 1))\
        * tcount_openbeamexp_noisy / np.sum(mean_openbeamexp)
        #* tcount_openbeamexp_noisy / np.sum(mean_openbeamexp) * 0.78
    op = ((np.zeros_like(data)+1)*mean_openbeamexp).transpose((2, 0, 1))\
        * tcount_openbeamexp_noisy/np.sum(mean_openbeamexp)
    print("FUCK", sample.shape, op.shape, data.shape)
    plt.plot((sample/op)[:, 40, 24])
    plt.show()
    x = np.arange(_shape[0])*40.+7240
    #output = np.vstack((x, np.random.poisson(sample[:, 40, 24]/10.)/np.random.poisson(op[:, 40, 24]/10.)))
    output = np.vstack((x, sample[:, 40, 24]/op[:, 40, 24]))
    #output = np.vstack((x, data[40, 24, :]))
    output = np.vstack((output, np.zeros_like(x) + 0.01))
    np.savetxt('check_simuspectrum.txt', output.T)

    # for median filter
    k = 24
    w = sampleexp.shape[0]
    idx = np.fromfunction(lambda i, j: i + j, (k, w), dtype=np.int) - k // 2
    idx[idx < 0] = 0
    idx[idx > w - 1] = w - 1
    #
    smoothed = np.median(sampleexp[idx, 40-8, 24-16], axis=0)
    fig, ax = plt.subplots(2, 2)
    for didx, (spectrum, label) in enumerate(zip([sampleexp_noisy[:, 40-8, 24-16],
                                             np.random.poisson(
                                                 sample[:, 24-16, 40-8])],
                                             ['sampleexp_noisy', 'samplesim'])):
    #for didx, (spectrum, label) in enumerate(zip([openbeamexp_noisy[:, 40, 24], np.random.poisson(op[:, 24, 40])],
    #                                             ['openbeamexp_noisy', 'op'])):
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
    plt.savefig('noise_levels_sample_exp_simu.png')
    plt.show()

    fig, ax = plt.subplots(4, 2, figsize=(10, 15))
    for didx, (spectrum, label) in enumerate(zip([sampleexp[:, 40-8, 24-16],
                                                  sampleexp_noisy[:, 40-8, 24-16],
                                                  openbeamexp[:, 40-8, 24-16],
                                                  openbeamexp_noisy[:, 40-8, 24-16]
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
    plt.savefig('noise_levels_sample_exp_full_part.png')
    plt.show()

    fig, ax = plt.subplots(4, 2, figsize=(10, 10))
    _op = np.random.poisson(op[:, :56-32, :])
    _openbeamexp_noisy = openbeamexp_noisy.transpose((0, 2, 1))
    vmin1 = np.min(np.vstack((_op[:, 4, :], _openbeamexp_noisy[:, 4, :])))
    vmax1 = np.max(np.vstack((_op[:, 4, :], _openbeamexp_noisy[:, 4, :])))
    vmin2 = np.min(np.vstack((_op[:, :, 32], _openbeamexp_noisy[:, :, 32])))
    vmax2 = np.max(np.vstack((_op[:, :, 32], _openbeamexp_noisy[:, :, 32])))
    vmin3 = np.min(np.vstack((_op.sum(axis=0), _openbeamexp_noisy.sum(axis=0))))
    vmax3 = np.max(np.vstack((_op.sum(axis=0), _openbeamexp_noisy.sum(axis=0))))
    ymax = np.max(np.vstack((_op[:, 4, 32], _openbeamexp_noisy[:, 4, 32])))
    ymin = np.min(np.vstack((_op[:, 4, 32], _openbeamexp_noisy[:, 4, 32])))
    for didx, (image, name) in enumerate(zip([_op, _openbeamexp_noisy],
                                             ['simulated_openbeam',
                                                 'expt_openbeam_noisy'])):
        ax[0, didx].set_title(name + ' crosssection at y=4 ch')
        ax[0, didx].set_xlabel(r'$\lambda$ / ch')
        ax[0, didx].set_ylabel('y / ch')
        ax[0, didx].imshow(image[:, 20-16, :].T, vmin=vmin1, vmax=vmax1, aspect=5)
        ax[1, didx].set_title(name + ' crosssection at x=32 ch')
        ax[1, didx].set_xlabel(r'$\lambda$ / ch')
        ax[1, didx].imshow(image[:, :, 40-8].T, vmin=vmin2, vmax=vmax2, aspect=10)
        ax[1, didx].set_ylabel('x / ch')
        ax[2, didx].set_title(name + ' projected onto xy')
        ax[2, didx].imshow(image.sum(axis=0), vmin=vmin3, vmax=vmax3)
        ax[2, didx].set_xlabel('x / ch')
        ax[2, didx].set_ylabel('y / ch')
        ax[2, didx].axhline(4, color='r', linestyle='--', lw=1)
        ax[2, didx].axvline(32, color='r', linestyle='--', lw=1)
        ax[3, didx].set_title(name + ' intensity at y=4 ch, x=32 ch')
        ax[3, didx].plot(image[:, 20-16, 40-8])
        ax[3, didx].set_ylim([ymin, ymax])
        ax[3, didx].set_xlabel(r'$\lambda$ / ch')
        ax[3, didx].set_ylabel('Neutron Count')
    plt.tight_layout()
    plt.savefig('openbeam_data_expt_simu.png')
    plt.show()
    fig, ax = plt.subplots(4, 2, figsize=(10, 10))
    _sample = np.random.poisson(sample[:, :56-32, :])
    _sampleexp_noisy = sampleexp_noisy.transpose((0, 2, 1))
    vmin1 = np.min(np.vstack((_sample[:, 4, :], _sampleexp_noisy[:, 4, :])))
    vmax1 = np.max(np.vstack((_sample[:, 4, :], _sampleexp_noisy[:, 4, :])))
    vmin2 = np.min(np.vstack((_sample[:, :, 32], _sampleexp_noisy[:, :, 32])))
    vmax2 = np.max(np.vstack((_sample[:, :, 32], _sampleexp_noisy[:, :, 32])))
    vmin3 = np.min(np.vstack((_sample.sum(axis=0), _sampleexp_noisy.sum(axis=0))))
    vmax3 = np.max(np.vstack((_sample.sum(axis=0), _sampleexp_noisy.sum(axis=0))))
    ymax = np.max(np.vstack((_sample[:, 4, 32], _sampleexp_noisy[:, 4, 32])))
    ymin = np.min(np.vstack((_sample[:, 4, 32], _sampleexp_noisy[:, 4, 32])))
    for didx, (image, name) in enumerate(zip([_sample, _sampleexp_noisy],
                                             ['simulated_sample',
                                                 'expt_sample_noisy'])):
        ax[0, didx].set_title(name + ' crosssection at y=4 ch')
        ax[0, didx].set_xlabel(r'$\lambda$ / ch')
        ax[0, didx].set_ylabel('y / ch')
        ax[0, didx].imshow(image[:, 20-16, :].T, vmin=vmin1, vmax=vmax1, aspect=5)
        ax[1, didx].set_title(name + ' crosssection at x=32 ch')
        ax[1, didx].set_xlabel(r'$\lambda$ / ch')
        ax[1, didx].imshow(image[:, :, 40-8].T, vmin=vmin2, vmax=vmax2, aspect=10)
        ax[1, didx].set_ylabel('x / ch')
        ax[2, didx].set_title(name + ' projected onto xy')
        ax[2, didx].imshow(image.sum(axis=0), vmin=vmin3, vmax=vmax3)
        ax[2, didx].set_xlabel('x / ch')
        ax[2, didx].set_ylabel('y / ch')
        ax[2, didx].axhline(4, color='r', linestyle='--', lw=1)
        ax[2, didx].axvline(32, color='r', linestyle='--', lw=1)
        ax[3, didx].set_title(name + ' intensity at y=4 ch, x=32 ch')
        ax[3, didx].plot(image[:, 20-16, 40-8])
        ax[3, didx].set_ylim([ymin, ymax])
        ax[3, didx].set_xlabel(r'$\lambda$ / ch')
        ax[3, didx].set_ylabel('Neutron Count')
    plt.tight_layout()
    plt.savefig('sample_data_expt_simu.png')
    plt.show()
    #plt.plot(op[300, 20-16, :], marker='o')
    #plt.show()
    #plt.plot(op[300, :56-32, 50-8], marker='o')
    #plt.show()
    sample = sample.transpose((1, 2, 0))
    op = op[:, :56, :].transpose((1, 2, 0))
    #plt.hist(sample.sum(axis=1).sum(axis=1))
    #plt.show()
    sim_datasets = {}
    sim_datasets['sample'] = sample
    sim_datasets['op'] = op
    datasets['x'] = x
    #with open('/home/kazu/desktop/240424/bi2d/bi2d_testbcc_simudata_rev2_lim.pkl', 'wb') as f:
    #    pickle.dump(sim_datasets, f, 4)


def gather_bi2d_only_cond(timescale=600, nidx=100):
    for iniidx in range(nidx):
        with open('/home/kazu/desktop/240424/bi2d/' + str(iniidx) + '_rev.pkl',
                  'rb') as f:
            data = pickle.load(f)
        cond = np.isnan(data).sum(axis=1).sum(axis=1) == 0
        if iniidx == 0:
            condt = cond
        else:
            condt = np.concatenate([condt, cond])
    param_sets_sets = load_param_sets_sets(param_sets_sets_file='param_sets_sets_bccrev.pkl')
    param_sets_sets_disel = []
    for pssidx in range(30000):
        if not condt[pssidx]:
            param_sets_sets_disel.append(param_sets_sets[pssidx])

    with open('/home/kazu/desktop/240424/param_sets_sets_bccrev_disel.pkl', 'wb') as f:
        pickle.dump(param_sets_sets_disel, f, 4)


def gather_bi2d(timescale=600, nidx=100):
    #with open('/home/kazu/desktop/240424/uNID_data_KO/433/openbeam.pkl',
    #          'rb') as f:
    #    openbeamexp = pickle.load(f)
    #    openbeamexp_noisy = pickle.load(f)
    #    etc = pickle.load(f)
    #mean_openbeamexp = openbeamexp.sum(axis=1).sum(axis=1)
    #tcount_openbeamexp_noisy = openbeamexp_noisy.sum(axis=0)
    for iniidx in range(nidx):
        print(iniidx)
        with open('/home/kazu/desktop/240424/bi2d/' + str(iniidx) + '_rev3.pkl',
                  'rb') as f:
            data = pickle.load(f)
        #print(data.shape)
        #plt.imshow(data[0, :, :])
        #plt.plot(data[0, 10, :])
        #plt.show()
        cond = np.isnan(data).sum(axis=1).sum(axis=1) == 0
        data = data[cond]
        if iniidx == 0:
            datat = data
            condt = cond
        else:
            datat = np.vstack((datat, data))
            condt = np.concatenate([condt, cond])
    #datat_noisy = get_noisydata(datat, x, timescale)
    datasets = {}
    datasets['target'] = datat*timescale
    #datasets['noisy'] = datat_noisy
    #datasets['x'] = x
    with open('/home/kazu/desktop/240424/bi2d/bi2d_testbcc_rev3.pkl', 'wb') as f:
        pickle.dump(datasets, f, 4)
        #pickle.dump(condt, f, 4)


def gather_bi3d(timescale=600, nidx=100):
    datat = []
    for iniidx in range(nidx):
        print(iniidx)
        with open('/home/kazu/desktop/240424/bi3d/bi2d_test_.pkl.'
                  + str(iniidx), 'rb') as f:
            data = pickle.load(f)
            x = pickle.load(f)
        data = data[np.isnan(data).sum(axis=1).sum(axis=1).sum(axis=1) == 0]
        if iniidx == 0:
            datat = data
        else:
            datat = np.vstack((datat, data))
    datat_noisy = get_noisydata(datat, x, timescale)
    datasets = {}
    datasets['target'] = datat*timescale
    datasets['noisy'] = datat_noisy
    datasets['x'] = x
    print('start to write data into a file')
    with open('/home/kazu/desktop/240424/bi3d/bi3d_test.pkl', 'wb') as f:
        pickle.dump(datasets, f, 4)


def select_bi2d(timescale=600):
    with open('/home/kazu/desktop/240424/bi2d/bi2d_test.pkl', 'rb') as f:
        data = pickle.load(f)
        x = pickle.load(f)
    print(data.shape)
    datat = data[np.isnan(data).sum(axis=1).sum(axis=1) == 0]
    datat_noisy = get_noisydata(datat, x, timescale)
    datasets = {}
    datasets['target'] = datat*timescale
    datasets['noisy'] = datat_noisy
    datasets['x'] = x
    with open('/home/kazu/desktop/240424/bi2d/bi2d_sel.pkl', 'wb') as f:
        pickle.dump(datasets, f, 4)


def check_data():
    with open('/home/kazu/desktop/240424/bi2d/bi2d_test.pkl', 'rb') as f:
        datasets = pickle.load(f)
    print(datasets['target'].shape)
    for i in range(5):
        _orig = datasets['target'][i]
        _noise = datasets['noisy'][i]
        bi2d.plot(_noise, _orig, _orig, savefig='test_bi_'+str(i)+'.png')


#gather_bi2d(timescale=50)
#cycles(ns=40000, dim=1)
#cycles_mpi(ns=240, dim=2, pklfile='param_sets_sets_2dlarge.pkl')
#divide_paramdata()
#mpi_parallel_computation(pklfile='param_sets_sets_2dlarge.pkl')
#single_computation()
#gather_bi2d(timescale=1, nidx=400)
#gather_bi2d_only_cond(timescale=1, nidx=300)
synthesize_bi2ddata()
#gather_bi3d(timescale=50, nidx=125)
#select_bi2d()
#check_data()
#crude_parallel_computation()
