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
from rits_fit_kt import get_sim_spectrum
sys.path.append("/home/kazu/denoise")
import bi2d
np.random.seed(129)


def get_params(rpkl=None, wpkl=None):
    if rpkl:
        with open(rpkl, 'rb') as f:
            param_sets = pickle.load(f)
    else:
        param_sets = {}
        param_sets['param_name'] = ["lattice constance"]
        param_sets['string'] = ['LATCON']
        param_sets['lims_mean'] = [[2.86, 2.87]]
        param_sets['lims_scale'] = [[0.0, 0.00005]]
        param_sets['lims_xlim'] = [[5., 25.]]
        param_sets['param_name'].append("crsytallite size")
        param_sets['string'].append("CRSIZE")
        param_sets['lims_mean'].append([0.186038, 1.31959])
        param_sets['lims_scale'].append([0.0, 0.01])
        param_sets['lims_xlim'].append([10., 30.])
        param_sets['param_name'].append("microstrain in Jorgensen func")
        param_sets['string'].append("MICRST")
        param_sets['lims_mean'].append([1.0, 15.7124])
        param_sets['lims_scale'].append([0.1, 2.0])
        param_sets['lims_xlim'].append([10., 30.])
        param_sets['param_name'].append("March-Dollase Coefficient")
        param_sets['string'].append("MDCOEF")
        param_sets['lims_mean'].append([1.33567, 2.53587])
        param_sets['lims_scale'].append([0.01, 1.0])
        param_sets['lims_xlim'].append([10., 30.])
        param_sets['param_name'].append("Projected Atomic Number Density")
        param_sets['string'].append("PRODEN")
        param_sets['lims_mean'].append([1.33478, 3.96254])
        param_sets['lims_scale'].append([0.05, 2.0])
        param_sets['lims_xlim'].append([6., 10.])
        if wpkl:
            with open(wpkl, 'wb') as f:
                pickle.dump(param_sets, f, 4)
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
                                  xlim=xlim, xsize=96)
        elif dim == 2:
            _params = draw_sample2d(mean=mean, numsample=1, scale=scale,
                                    xlim=xlim, xsize=72, rsize=24)
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
                                    scale=scale, xlim=xlim, xsize=72,
                                    ysize=24, rsize=12)
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
    param_sets = get_params(wpkl='param_sets_bccrev2.pkl')
    for ins in range(ns):
        param_sets_sets.append(draw_params(param_sets, dim=dim))
    if dim == 1:
        with open('param_sets_sets_bccrev5.pkl', 'wb') as f:
            pickle.dump(param_sets_sets, f, 4)
    elif dim == 2:
        with open('param_sets_sets_bccrev2_2d.pkl', 'wb') as f:
            pickle.dump(param_sets_sets, f, 4)


def cycles_mpi(ns=10, dim=2, pklfile='param_sets_sets_bccrev2_2d.pkl'):
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
    param_sets = get_params(rpkl='param_sets_bccrev2.pkl')
    for ins in range(rank*(ns//psize), (rank+1)*(ns//psize)):
        print('rank=', rank, 'ins=', ins)
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
        print("run_rits pidx", pidx,'/',params[0].size)
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
        #bi2d_true = bi2d_true.reshape((params[0].shape[0], params[0].shape[1], x.shape[0]))
        bi2d_true = bi2d_true.reshape((params[0].shape[1], params[0].shape[0], x.shape[0]))
    return bi2d_true, x


def get_noisydata(bi2d, x, timescale):
    return np.random.poisson(bi2d*timescale*x)/x


# routine for crude parallel computations of rits
def crude_parallel_computation(pklfile='param_sets_sets_bccrev5.pkl'):
    inino = int(sys.argv[1])
    inpfile = 'rits_initial.inp.' + str(inino)
    param_sets_sets = load_param_sets_sets(param_sets_sets_file=pklfile)
    for pid, param_sets in enumerate(param_sets_sets[inino*100:inino*100+100]):
        bi2d_true, x = run_rits(param_sets['params'],  param_sets['string'],
                                inpfile=inpfile)
        if pid == 0:
            bi2dt = np.zeros((100, bi2d_true.shape[0], bi2d_true.shape[1]))
        bi2dt[pid] = bi2d_true
    with open('/home/kazu/desktop/240424/bi2d/' + str(inino) + '_rev5large.pkl', 'wb') as f:
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
    with open('bi2dsingle.pkl.'+str(iniidx), 'wb') as f:
        pickle.dump(bi3d_true, f, 4)
        pickle.dump(x, f, 4)


def synthesize_bi2ddata():
    #with open('/home/kazu/desktop/240424/suuNID/433/openbeamstrd_openbaemnoisylowcount.pkl',
    with open('/home/kazu/desktop/240424/suuNID/333/openbeamstrd_openbaemnoisylowcount.pkl',
              'rb') as f:
        # 3020	#4#	---23780	#5#	---
        #initof = (3020-20)//40 
        #fintof = (23780-20)//40 + 1
        #openbeamexp = pickle.load(f)[initof:fintof, 2:34, 2:22]
        #openbeamexp_noisy = pickle.load(f)[initof:fintof, 2:34, 2:22]
        #sampleexp = pickle.load(f)[initof:fintof, 2:34, 2:22]
        #sampleexp_noisy = pickle.load(f)[initof:fintof, 2:34, 2:22]
        openbeamexp = pickle.load(f)[:, 2:98, 2:62]
        openbeamexp_noisy = pickle.load(f)[:, 2:98, 2:62]
        sampleexp = pickle.load(f)[:, 2:98, 2:62]
        sampleexp_noisy = pickle.load(f)[:, 2:98, 2:62]
        print("CHK", openbeamexp.shape)
    with open('/home/kazu/desktop/240424/bi2d/bi2d_testbcc_rev5large.pkl',
              'rb') as f:
        datasets = pickle.load(f)
    data = datasets['target'][:, :, :]
    print(data.shape)
    _shape = openbeamexp.shape
    tcount_openbeamexp_noisy = np.tile(openbeamexp_noisy.sum(axis=0),
                                       (1, 3440)).T[:data.shape[0]]
    mean_openbeamexp = openbeamexp.transpose((2, 1, 0))\
        .sum(axis=0).sum(axis=0)/_shape[1]/_shape[2]*1.15
    plt.plot(mean_openbeamexp)
    plt.show()
    # correction due to the different # of kickers is added.
    sample = (data*mean_openbeamexp).transpose((2, 0, 1))\
        * tcount_openbeamexp_noisy / np.sum(mean_openbeamexp)*199981/113494
    op = ((np.zeros_like(data)+1)*mean_openbeamexp).transpose((2, 0, 1))\
        * tcount_openbeamexp_noisy/np.sum(mean_openbeamexp)
    print("FUCK", sample.shape, op.shape, data.shape)
    vpos = 24*3
    hpos = 20
    plt.plot((sample/op)[:, hpos, vpos])
    plt.show()
    x = np.arange(_shape[0])*40.+3020
    output = np.vstack((x, sample[:, hpos, vpos]/op[:, hpos, vpos]))
    output = np.vstack((output, np.zeros_like(x) + 0.01))
    np.savetxt('check_simuspectrum.txt', output.T)

    # for median filter
    k = 24
    w = sampleexp.shape[0]
    idx = np.fromfunction(lambda i, j: i + j, (k, w), dtype=np.int) - k // 2
    idx[idx < 0] = 0
    idx[idx > w - 1] = w - 1
    print(sampleexp.shape)
    smoothed = np.median(sampleexp[idx, vpos, hpos], axis=0)
    fig, ax = plt.subplots(2, 2)
    for didx, (sp, lab) in enumerate(zip([sampleexp_noisy[:, vpos, hpos],
                                     np.random.poisson(sample[:, hpos, vpos])],
                                     ['sampleexp_noisy', 'samplesim'])):
    #for didx, (sp, lab) in enumerate(zip([openbeamexp_noisy[:, 40, 24], np.random.poisson(op[:, 24, 40])],
    #                                             ['openbeamexp_noisy', 'op'])):
        smoothed = np.median(sp[idx], axis=0)
        ax[0, didx].plot(sp, label=lab)
        ax[0, didx].plot(smoothed, label='medianfiltered')
        ax[0, didx].set_xlabel(r'$\lambda$ / ch')
        ax[0, didx].set_xlabel('TOF / ch')
        ax[0, didx].set_ylabel('Neutron Count')
        ax[0, didx].legend()
        ax[1, didx].plot(sp - smoothed)
        ax[1, didx].set_ylabel(lab + '-medianfiltered')
        ax[1, didx].set_xlabel('TOF / ch')
    plt.tight_layout()
    #plt.savefig('noise_levels_sample_exp_simulong.png')
    plt.show()

    fig, ax = plt.subplots(4, 2, figsize=(10, 15))
    for didx, (spectrum, label) in enumerate(zip([sampleexp[:, vpos, hpos],
                                             sampleexp_noisy[:, vpos, hpos],
                                             openbeamexp[:, vpos, hpos],
                                             openbeamexp_noisy[:, vpos, hpos]],
                                             ['sample', 'sample_part',
                                              'openbeam', 'openbeam_part'])):
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
    plt.savefig('noise_levels_sample_exp_full_partlong.png')
    plt.show()

    fig, ax = plt.subplots(4, 2, figsize=(10, 10))
    _op = np.random.poisson(op[:, :_shape[2], :])
    _openbeamexp_noisy = openbeamexp_noisy.transpose((0, 2, 1))
    vmin1 = np.min(np.vstack((_op[:, hpos, :],
                              _openbeamexp_noisy[:, hpos, :])))
    vmax1 = np.max(np.vstack((_op[:, hpos, :],
                              _openbeamexp_noisy[:, hpos, :])))
    vmin2 = np.min(np.vstack((_op[:, :, vpos],
                              _openbeamexp_noisy[:, :, vpos])))
    vmax2 = np.max(np.vstack((_op[:, :, vpos],
                              _openbeamexp_noisy[:, :, vpos])))
    vmin3 = np.min(np.vstack((_op.sum(axis=0),
                              _openbeamexp_noisy.sum(axis=0))))
    vmax3 = np.max(np.vstack((_op.sum(axis=0),
                              _openbeamexp_noisy.sum(axis=0))))
    ymax = np.max(np.vstack((_op[:, hpos, vpos],
                             _openbeamexp_noisy[:, hpos, vpos])))
    ymin = np.min(np.vstack((_op[:, hpos, vpos],
                             _openbeamexp_noisy[:, hpos, vpos])))
    for didx, (image, name) in enumerate(zip([_op, _openbeamexp_noisy],
                                         ['simulated_openbeam',
                                          'expt_openbeam_noisy'])):
        ax[0, didx].set_title(f'{name} crosssection at y={hpos} ch')
        ax[0, didx].set_xlabel(r'$\lambda$ / ch')
        ax[0, didx].set_ylabel('y / ch')
        ax[0, didx].imshow(image[:, hpos, :].T, vmin=vmin1, vmax=vmax1,
                           aspect=5)
        ax[1, didx].set_title(f'{name} crosssection at x={vpos} ch')
        ax[1, didx].set_xlabel(r'$\lambda$ / ch')
        ax[1, didx].imshow(image[:, :, vpos].T, vmin=vmin2, vmax=vmax2,
                           aspect=10)
        ax[1, didx].set_ylabel('x / ch')
        ax[2, didx].set_title(f'{name} projected onto xy')
        ax[2, didx].imshow(image.sum(axis=0), vmin=vmin3, vmax=vmax3)
        ax[2, didx].set_xlabel('x / ch')
        ax[2, didx].set_ylabel('y / ch')
        ax[2, didx].axhline(hpos, color='r', linestyle='--', lw=1)
        ax[2, didx].axvline(vpos, color='r', linestyle='--', lw=1)
        ax[3, didx].set_title(f'{name} intensity at y={hpos} ch, x={vpos} ch')
        ax[3, didx].plot(image[:, hpos, vpos])
        ax[3, didx].set_ylim([ymin, ymax])
        ax[3, didx].set_xlabel(r'$\lambda$ / ch')
        ax[3, didx].set_ylabel('Neutron Count')
    plt.tight_layout()
    #plt.savefig('openbeam_data_expt_simulong.png')
    plt.show()
    fig, ax = plt.subplots(4, 2, figsize=(10, 10))
    _sample = np.random.poisson(sample[:, :_shape[2], :])
    _sampleexp_noisy = sampleexp_noisy.transpose((0, 2, 1))
    vmin1 = np.min(np.vstack((_sample[:, hpos, :],
                   _sampleexp_noisy[:, hpos, :])))
    vmax1 = np.max(np.vstack((_sample[:, hpos, :],
                   _sampleexp_noisy[:, hpos, :])))
    vmin2 = np.min(np.vstack((_sample[:, :, vpos],
                   _sampleexp_noisy[:, :, vpos])))
    vmax2 = np.max(np.vstack((_sample[:, :, vpos],
                   _sampleexp_noisy[:, :, vpos])))
    vmin3 = np.min(np.vstack((_sample.sum(axis=0),
                   _sampleexp_noisy.sum(axis=0))))
    vmax3 = np.max(np.vstack((_sample.sum(axis=0),
                   _sampleexp_noisy.sum(axis=0))))
    ymax = np.max(np.vstack((_sample[:, hpos, vpos],
                  _sampleexp_noisy[:, hpos, vpos])))
    ymin = np.min(np.vstack((_sample[:, hpos, vpos],
                  _sampleexp_noisy[:, hpos, vpos])))
    for didx, (image, name) in enumerate(zip([_sample, _sampleexp_noisy],
                                             ['simulated_sample',
                                                 'expt_sample_noisy'])):
        ax[0, didx].set_title(f'{name} crosssection at y={hpos} ch')
        ax[0, didx].set_xlabel(r'$\lambda$ / ch')
        ax[0, didx].set_ylabel('y / ch')
        ax[0, didx].imshow(image[:, hpos, :].T, vmin=vmin1, vmax=vmax1,
                           aspect=5)
        ax[1, didx].set_title(f'{name} crosssection at x={vpos} ch')
        ax[1, didx].set_xlabel(r'$\lambda$ / ch')
        ax[1, didx].imshow(image[:, :, vpos].T, vmin=vmin2, vmax=vmax2,
                           aspect=10)
        ax[1, didx].set_ylabel('x / ch')
        ax[2, didx].set_title(f'{name} projected onto xy')
        ax[2, didx].imshow(image.sum(axis=0), vmin=vmin3, vmax=vmax3)
        ax[2, didx].set_xlabel('x / ch')
        ax[2, didx].set_ylabel('y / ch')
        ax[2, didx].axhline(hpos, color='r', linestyle='--', lw=1)
        ax[2, didx].axvline(vpos, color='r', linestyle='--', lw=1)
        ax[3, didx].set_title(f'{name} intensity at y={hpos} ch, x={vpos} ch')
        ax[3, didx].plot(image[:, hpos, vpos])
        ax[3, didx].set_ylim([ymin, ymax])
        ax[3, didx].set_xlabel(r'$\lambda$ / ch')
        ax[3, didx].set_ylabel('Neutron Count')
    plt.tight_layout()
    #plt.savefig('sample_data_expt_simulong.png')
    plt.show()
    #plt.plot(op[300, 20-16, :], marker='o')
    #plt.show()
    #plt.plot(op[300, :56-32, 50-8], marker='o')
    #plt.show()
    sample = sample.transpose((1, 2, 0))
    op = op[:, :_shape[2], :].transpose((1, 2, 0))
    #plt.hist(sample.sum(axis=1).sum(axis=1))
    #plt.show()
    sim_datasets = {}
    sim_datasets['sample'] = sample
    sim_datasets['op'] = op
    datasets['x'] = x
    with open('/home/kazu/desktop/240424/bi2d/' +
              'bi2d_testbcc_simudata_rev5large_lim.pkl', 'wb') as f:
        pickle.dump(sim_datasets, f, 4)


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
        with open('/home/kazu/desktop/240424/bi2d/' + str(iniidx) + '_rev5long.pkl',
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
            #condt = cond
        else:
            datat = np.vstack((datat, data))
            #condt = np.concatenate([condt, cond])
    #datat_noisy = get_noisydata(datat, x, timescale)
    datasets = {}
    datasets['target'] = datat*timescale
    #datasets['noisy'] = datat_noisy
    #datasets['x'] = x
    with open('/home/kazu/desktop/240424/bi2d/bi2d_testbcc_rev5long.pkl', 'wb') as f:
        pickle.dump(datasets, f, 4)
        #pickle.dump(condt, f, 4)


def gather_bi2d_mpi(timescale=600, nidx=100):
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    psize = comm.Get_size()
    if rank == 0:
        print(f'psize:{psize}')
    if nidx % psize == 0:
        if rank == 0:
            print("divided")
        numregular = nidx // psize
        ininum = rank*numregular
        finnum = (rank+1)*numregular
    else:
        if rank == 0:
            print("not divided")
        if rank <= (nidx % psize) - 1:
            numregular = nidx // psize + 1
            ininum = rank*numregular
            finnum = (rank+1)*numregular
        else:
            numregular = nidx // psize + 1
            numregular2 = (nidx - (nidx % psize)) // psize
            ininum = (nidx % psize) * numregular + (rank - (nidx % psize)
                                                    )*numregular2
            finnum = (nidx % psize) * numregular + (rank - (nidx % psize) + 1
                                                    )*numregular2
    print(f'ininum:{ininum}, finnum:{finnum}, rank:{rank}')
    for iniidx in range(ininum, finnum):
        if iniidx % 10 == 0:
            print(f'iniidx {iniidx} @ rank {rank}')
        with open(f'/home/kazu/desktop/240424/bi2d/{iniidx}_rev5large.pkl',
                  'rb') as f:
            data = pickle.load(f)
        cond = np.isnan(data).sum(axis=1).sum(axis=1) == 0
        data = data[cond]
        if iniidx == ininum:
            datat = data
        else:
            datat = np.vstack((datat, data))
    with open('/home/kazu/desktop/240424/bi2d/' +
              f'bi2d_testbcc_rev5large_{rank}.pkl', 'wb') as f:
        pickle.dump(datat, f, 4)
    comm.barrier()
    if rank == 0:
        for ridx in range(psize):
            with open('/home/kazu/desktop/240424/bi2d/' +
                      f'bi2d_testbcc_rev5large_{ridx}.pkl', 'rb') as f:
                data = pickle.load(f)
            if ridx == 0:
                datat = data
            else:
                datat = np.vstack((datat, data))
        datasets = {}
        datasets['target'] = datat*timescale
        with open('/home/kazu/desktop/240424/bi2d/'
                  + 'bi2d_testbcc_rev5large.pkl', 'wb') as f:
            pickle.dump(datasets, f, 4)


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
#cycles(ns=2, dim=2)
#cycles_mpi(ns=1200, dim=2, pklfile='param_sets_sets_bccrev2_2d.pkl')
#divide_paramdata()
#mpi_parallel_computation(pklfile='param_sets_sets_2dlarge.pkl')
single_computation(pklfile='param_sets_sets_bccrev2_2d.pkl')
#gather_bi2d(timescale=1, nidx=300)
#gather_bi2d_mpi(timescale=1, nidx=300)
#gather_bi2d_only_cond(timescale=1, nidx=300)
#synthesize_bi2ddata()
#gather_bi3d(timescale=50, nidx=125)
#select_bi2d()
#check_data()
#crude_parallel_computation()
