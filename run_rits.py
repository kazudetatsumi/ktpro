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
from rits_fit_kt import get_sim_spectrum, get_sim_edgespectrum
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
        param_sets['lims_xlim'] = [[25., 45.]]
        param_sets['param_name'].append("crsytallite size")
        param_sets['string'].append("CRSIZE")
        param_sets['lims_mean'].append([0.186038, 1.31959])
        param_sets['lims_scale'].append([0.0, 0.01])
        param_sets['lims_xlim'].append([30., 50.])
        param_sets['param_name'].append("microstrain in Jorgensen func")
        param_sets['string'].append("MICRST")
        param_sets['lims_mean'].append([1.0, 15.7124])
        param_sets['lims_scale'].append([0.1, 2.0])
        param_sets['lims_xlim'].append([30., 50.])
        param_sets['param_name'].append("March-Dollase Coefficient")
        param_sets['string'].append("MDCOEF")
        param_sets['lims_mean'].append([1.33567, 2.53587])
        param_sets['lims_scale'].append([0.01, 1.0])
        param_sets['lims_xlim'].append([30., 50.])
        param_sets['param_name'].append("Projected Atomic Number Density")
        param_sets['string'].append("PRODEN")
        param_sets['lims_mean'].append([1.33478, 3.96254])
        param_sets['lims_scale'].append([0.05, 2.0])
        param_sets['lims_xlim'].append([10., 14.])
        if wpkl:
            with open(wpkl, 'wb') as f:
                pickle.dump(param_sets, f, 4)
    return param_sets


def get_params_edge(rpkl=None, wpkl=None):
    if rpkl:
        with open(rpkl, 'rb') as f:
            param_sets = pickle.load(f)
    else:
        param_sets = {}
        param_sets['param_name'] = ["Tend"]
        param_sets['string'] = ['Te']
        param_sets['lims_mean'] = [[0.7, 1.0]]
        param_sets['lims_scale'] = [[0.001, 0.03]]
        param_sets['lims_xlim'] = [[10., 30.]]
        param_sets['param_name'].append("Alpha")
        param_sets['string'].append("alpha")
        param_sets['lims_mean'].append([0.9875, 1.0125])
        param_sets['lims_scale'].append([0.001, 0.0025])
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
        #param_sets['lims_mean'].append([2.025, 2.038])
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
                                    xlim=xlim, xsize=96, rsize=56)
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
        while True:
            mean = rg.uniform(low=lims_mean[0], high=lims_mean[1])
            scale = rg.uniform(low=lims_scale[0], high=lims_scale[1])
            xlim = rg.uniform(low=lims_xlim[0], high=lims_xlim[1])
            if dim == 1:
                _params = draw_sample(mean=mean, numsample=1, scale=scale,
                                      xlim=xlim)
            elif dim == 2:
                _params, _kerneltype = draw_sample2d_mpi(
                        rg, mean=mean, numsample=1, scale=scale, xlim=xlim,
                        xsize=72, ysize=192, rsize=1)
            if pidx != 3 and np.min(_params) > 0.:
                break
            if pidx == 3 and np.min(_params) > 0.35:
                break
        #if pidx==0:
        #    print(param_sets['param_name'][pidx], param_sets['lims_mean'][pidx], mean, _params)
        #if np.min(_params) < 0.:
        #    _params -= np.min(_params)
        if pidx == 0:
            _param_sets = copy.deepcopy(param_sets)
            _param_sets['mean'] = [mean]
            _param_sets['scale'] = [scale]
            _param_sets['xlim'] = [xlim]
            _param_sets['params'] = [_params[0]]
        else:
            _param_sets['mean'].append(mean)
            _param_sets['scale'].append(scale)
            _param_sets['xlim'].append(xlim)
            _param_sets['params'].append(_params[0])
    return _param_sets


def draw_params_edge_mpi(rg, param_sets, dim=1, kerneltype='square'):
    if kerneltype == 'random':
        kerneltypes = ['square', 'single', 'third', 'fifth']
        _kerneltype = kerneltypes[rg.integers(4)]
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


def cycles(ns=2, dim=1):
    param_sets_sets = []
    param_sets = get_params(wpkl='param_sets_bccrev5.pkl')
    for ins in range(ns):
        param_sets_sets.append(draw_params(param_sets, dim=dim))
    if dim == 1:
        with open('param_sets_sets_bccrev5.pkl', 'wb') as f:
            pickle.dump(param_sets_sets, f, 4)
    elif dim == 2:
        with open('param_sets_sets_bccrev5_2d.pkl', 'wb') as f:
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


def cycles_edge_mpi(ns=10, dim=2, kerneltype='square',
                    pklfile='param_sets_sets_bccrev2_2d_edge.pkl'):
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
    param_sets = get_params_edge(rpkl='param_sets_bccrev2_edge.pkl')
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


def cycles_mpi_div(nss=10, ns=10, dim=2,
                   orgpklfile='param_sets_sets_bccrev2_2d.pkl'):
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    psize = comm.Get_size()
    import gc
    # for restart job
    param_sets = get_params(rpkl='param_sets_bccrev2.pkl')
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
        #else:
        #    from numpy.random import Generator, PCG64, SeedSequence
        #    sg = SeedSequence(1234)
        #    ss = sg.spawn(psize)
        #    rg = Generator(PCG64(ss[rank]))
        _param_sets_sets = []
        for ins in range(rank*(ns//psize), (rank+1)*(ns//psize)):
            print('rank=', rank, 'ins=', ins)
            _param_sets_sets.append(draw_params_mpi(rg, param_sets, dim=dim))
        comm.barrier()
        param_sets_sets = comm.gather(_param_sets_sets, root=0)
        rg_sets = comm.gather(rg, root=0)
        if rank == 0:
            __param_sets_sets = [__cont for _cont in param_sets_sets for
                                 __cont in _cont]
            if dim == 1:
                with open(outpklfile, 'wb') as f:
                    pickle.dump(param_sets_sets, f, 4)
            elif dim == 2:
                with open(outpklfile, 'wb') as f:
                    pickle.dump(__param_sets_sets, f, 4)
                    pickle.dump(rg_sets, f, 4)
            print(datetime.datetime.now(), outpklfile, ' is saved')
        comm.barrier()


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
        #print("run_rits pidx", pidx,'/',params[0].size)
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
    param_sets = load_param_sets_sets(param_sets_sets_file=pklfile)[iniidx]
    inpfile = 'rits_initial.inp.' + str(iniidx)
    bi3d_true, x = run_rits(param_sets['params'],  param_sets['string'],
                            inpfile=inpfile)
    with open('bi2dsingle.pkl.'+str(iniidx), 'wb') as f:
        pickle.dump(bi3d_true, f, 4)
        pickle.dump(x, f, 4)


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
        mask[didx] = gaussian_filter(mask[didx], sigma=np.random.uniform(low=maskparams[1],
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
    #print('chk', np.min(a0, axis=(1,2)) < 0.)
    #print('chk', np.min(b0, axis=(1,2))<0.)
    #print('chk', np.min(ahkl, axis=(1,2))<0.)
    #print('chk', np.min(bhkl, axis=(1,2)) < 0.)
    #for i in range(8):
    #    plt.imshow(a0[i])
    #    plt.show()
    #    plt.plot(a0[i, 10])
    #    plt.plot(a0[i, 20])
    #    plt.plot(a0[i, 30])
    #    plt.show()



#def single_computation_div(pklfile='param_sets_sets_bccrev.pkl'):
#    # USAGE: run_rits.py [iniidx, int] [tag, str]
#    # EXAMPLE: run_rits.py 0 3-1
#    iniidx = int(sys.argv[1])
#    tag = str(sys.argv[2])
#    pklfile = pklfile + "." + tag
#    param_sets = load_param_sets_sets(param_sets_sets_file=pklfile)[iniidx]
#    inpfile = 'rits_initial.inp.' + tag + "." + str(iniidx)
#    bi3d_true, x = run_rits(param_sets['params'],  param_sets['string'],
#                            inpfile=inpfile)
#    with open('bi2dsingle.pkl.' + tag + "." + str(iniidx), 'wb') as f:
#        pickle.dump(bi3d_true, f, 4)
#        pickle.dump(x, f, 4)


def single_edgecomputation_phantom(pklfile='paramimage_211_d2_5models.pkl'):
    # for constructing phantom data from the result of the rits edge fitting.
    with open(pklfile, 'rb') as f:
        paramimage = pickle.load(f)
    inpfile = 'edge_3.inp.phantom'
    bi3d_true, x = run_edge_phantom(paramimage, inpfile=inpfile)
    with open('bi2dsingle_211_d2_5models.pkl.phantom', 'wb') as f:
        pickle.dump(bi3d_true, f, 4)
        pickle.dump(x, f, 4)

def resample_mpi(testdata, pklfile, ns=10):
    # bootstrapping for a large number of counts by using mpi
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    psize = comm.Get_size()
    import gc

    from numpy.random import Generator, PCG64, SeedSequence
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
        sg = SeedSequence(1234)
        ss = sg.spawn(psize)
        rg = Generator(PCG64(ss[rank]))
    data = np.int32(testdata.flatten())
    datasize = data.shape[0]
    x = np.array([idx for idx in range(datasize)
                  for num_repeat in range(data[idx])], dtype=np.int32)
    rebuild_data_sets = []
    for ins in range(rank*(ns//psize), (rank+1)*(ns//psize)):
        Nb = rg.poisson(lam=len(x))
        idx = rg.integers(0, len(x), Nb, dtype=np.int32)
        xb = x[idx]
        _rebuild_data, bin_edge = np.histogram(xb, bins=datasize,
                                               range=(0, datasize))
        rebuild_data_sets.append(_rebuild_data)
    comm.barrier()
    rebuild_data_sets_sets = np.array(comm.gather(rebuild_data_sets, root=0))
    rg_sets = comm.gather(rg, root=0)
    if rank == 0:
        rebuild_data_sets_sets = rebuild_data_sets_sets.astype('int32')
        rebuild_data_sets_sets = rebuild_data_sets_sets.reshape((ns,) +
                                                                testdata.shape)
        if os.path.isfile(pklfile):
            with open(pklfile, 'rb') as f:
                _rebuild_data_sets_sets = pickle.load(f)
            rebuild_data_sets_sets = np.vstack((rebuild_data_sets_sets,
                                                _rebuild_data_sets_sets))
        with open(pklfile, 'wb') as f:
            pickle.dump(rebuild_data_sets_sets, f, 4)
            pickle.dump(rg_sets, f, 4)


def run_resample_mpi(ns=4):
    #with open('/home/kazu/desktop/240424/uNID_data_KO/211/openbeam.pkl',
    #          'rb') as f:
    #    openbeamexp = pickle.load(f)
    #    openbeamexp_noisy = pickle.load(f)
    #    sampleexp = pickle.load(f)
    #    sampleexp_noisy = pickle.load(f)
    with open('/home/kazu/desktop/240424/connect2d/sampled/bi3d_testbcc_simudata_rev2_lim_single_resize_full_255_d_phantom.pkl',
              'rb') as f:
        sample = pickle.load(f)
        op = pickle.load(f)
    resample_mpi(sample, 'sample_resample_d_phantom.pkl', ns=ns)
    resample_mpi(op, 'openbeam_resample_d_phantom.pkl', ns=ns)


def synthesize_bi3ddata():
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
    #datat = []
    nremainders = 0
    for iniidx in range(nidx):
        print(nremainders, iniidx)
        with open('bi2dsingle.pkl.' + str(iniidx), 'rb') as f:
            data = pickle.load(f)
            x = pickle.load(f)
            if np.sum(np.isnan(data)) == 0:
                nremainders += 1
        #data = data[np.isnan(data).sum(axis=1).sum(axis=1).sum(axis=1) == 0]
        #if iniidx == 0:
        #    datat = data
        #else:
        #    datat = np.vstack((datat, data))
    #datat_noisy = get_noisydata(datat, x, timescale)
    #datasets = {}
    #datasets['target'] = datat*timescale
    #datasets['noisy'] = datat_noisy
    #datasets['x'] = x
    #print('start to write data into a file')
    #with open('bi3d_test.pkl', 'wb') as f:
    #    pickle.dump(datasets, f, 4)


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
#cycles_mpi(ns=2560, dim=2, pklfile='param_sets_sets_bccrev2_2d_single_edge_MDCoeffrev.pkl')
#cycles_edge_mpi(ns=1280, dim=2, kerneltype='random', pklfile='param_sets_sets_bccrev2_2d_single_edge_true_edge_ktrand.pkl')
#cycles_mpi_div(nss=11, ns=2560, dim=2, orgpklfile='param_sets_sets_bccrev2_2d_single_edge_MDCoeffrev.pkl')
#cycles_edge_mpi_div(nss=11, ns=1280, dim=2, kerneltype='random', orgpklfile='param_sets_sets_bccrev2_2d_single_edge_true_edge_ktrand.pkl')
#divide_paramdata()
#mpi_parallel_computation(pklfile='param_sets_sets_2dlarge.pkl')
#single_computation(pklfile='param_sets_sets_bccrev2_2d_single_edge_MDCoeffrev.pkl')
#single_computation_div(pklfile='param_sets_sets_bccrev2_2d_single_edge_rsize.pkl')
#single_edgecomputation(pklfile='param_sets_sets_bccrev2_2d_single_edge_true_edge.pkl')
single_edgecomputation_phantom()
#gather_bi2d(timescale=1, nidx=300)
#gather_bi2d_mpi(timescale=1, nidx=300)
#gather_bi2d_only_cond(timescale=1, nidx=300)
#synthesize_bi3ddata()
#gather_bi3d(timescale=1, nidx=1200)
#select_bi2d()
#check_data()
#crude_parallel_computation()
#run_resample_mpi(ns=10)
