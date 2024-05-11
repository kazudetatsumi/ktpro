#!/usr/bin/env python
import numpy as np
import textfile
import os
import sys
import pickle
import matplotlib.pyplot as plt
import copy
sys.path.append("/home/kazu/ktpro")
from gp_nrca import draw_sample, draw_sample2d
#from gp_nrca_mpi import draw_sample2d
from rits_fit_kt import get_sim_spectrum
sys.path.append("/home/kazu/denoise")
import bi2d
np.random.seed(126)


def get_params():
    param_sets = {}
    param_sets['param_name'] = ["lattice constance"]
    param_sets['string'] = ['LATCON']
    param_sets['lims_mean'] = [[2.5, 4.6]]
    param_sets['lims_scale'] = [[0.0, 0.05]]
    param_sets['lims_xlim'] = [[10., 30.]]
    #45 Crsytallite size, 2.0〜9.0 ミクロン. エッジのジャンプ高さが低くなる。
    param_sets['param_name'].append("crsytallite size")
    param_sets['string'].append("CRSIZE")
    param_sets['lims_mean'].append([2.0, 9.0])
    param_sets['lims_scale'].append([0.1, 2.0])
    param_sets['lims_xlim'].append([10., 30.])
    #24 Crysallite size in Jorgensen function, エッジの傾きが変わる。
    # 0〜6ぐらい．
    param_sets['param_name'].append("crsytallite size in Jorgensen func")
    param_sets['string'].append("CRSIZJ")
    param_sets['lims_mean'].append([0.0, 6.0])
    param_sets['lims_scale'].append([0.1, 2.0])
    param_sets['lims_xlim'].append([10., 30.])
    #32 March-Dollase 係数,  1.5〜2.0 (蘇さんの歯車), 0.3〜0.9 (佐藤さんの溶接材)。
    # <1でビームに平行に成長、>1で垂直に成長。
    # 特定のhkl列のエッジジャンプが変わる．
    param_sets['param_name'].append("March-Dollase Coefficient")
    param_sets['string'].append("MDCOEF")
    param_sets['lims_mean'].append([0.3, 2.0])
    param_sets['lims_scale'].append([0.01, 0.1])
    param_sets['lims_xlim'].append([10., 30.])
    # 元素特性1
    param_sets['param_name'].append("Coherent Scattering Length")
    param_sets['string'].append("COHESL")
    param_sets['lims_mean'].append([9.4, 9.6])
    param_sets['lims_scale'].append([0.01, 0.1])
    param_sets['lims_xlim'].append([10., 30.])
    # 投影原子数密度
    param_sets['param_name'].append("Projected Atomic Number Density")
    param_sets['string'].append("PRODEN")
    param_sets['lims_mean'].append([6.3, 10.3])
    param_sets['lims_scale'].append([0.05, 0.5])
    param_sets['lims_xlim'].append([6., 10.])

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
                                    xlim=xlim, xsize=72, rsize=96)
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
        with open('param_sets_sets.pkl', 'wb') as f:
            pickle.dump(param_sets_sets, f, 4)
    elif dim == 2:
        with open('param_sets_sets_2d.pkl', 'wb') as f:
            pickle.dump(param_sets_sets, f, 4)


def load_param_sets_sets(param_sets_sets_file='param_sets_sets.pkl'):
    with open(param_sets_sets_file, 'rb') as f:
        param_sets_sets = pickle.load(f)
    return param_sets_sets


def run_rits(params, strings, inpfile='rits_initial.inp'):
    for pidx in range(len(params[0])):
        os.system('cp rits_initial.inp.temp '+inpfile)
        for sidx, string in enumerate(strings):
            _param = "{:.4f}".format(params[sidx][pidx])
            textfile.replace(inpfile, string, _param)
            if string == "COHESL" and pidx == 0:
                params2 = params[sidx]**2*4*np.pi/100.
            if string == "COHESL":
                _param2 = "{:.4f}".format(params2[pidx])
                textfile.replace(inpfile, "COHXSC", _param2)
        x, y = get_sim_spectrum(inpfile=inpfile)
        if pidx == 0:
            bi2d_true = np.zeros((params[0].shape[0], x.shape[0]))
        bi2d_true[pidx] = y
    return bi2d_true, x


def get_noisydata(bi2d, x, timescale):
    return np.random.poisson(bi2d*timescale*x)/x


#cycles(ns=10000)

# routine for crude parallel computations of rits
def crude_parallel_computation():
    inino = int(sys.argv[1])
    inpfile = 'rits_initial.inp.' + str(inino)
    param_sets_sets = load_param_sets_sets()
    for pid, param_sets in enumerate(param_sets_sets[inino*100:inino*100+100]):
        bi2d_true, x = run_rits(param_sets['params'],  param_sets['string'],
                                inpfile=inpfile)
        if pid == 0:
            bi2dt = np.zeros((100, bi2d_true.shape[0], bi2d_true.shape[1]))
        bi2dt[pid] = bi2d_true
    with open('/home/kazu/desktop/240424/bi2d/' + str(inino) + '.pkl', 'wb') as f:
        pickle.dump(bi2dt, f, 4)
        pickle.dump(x, f, 4)
    #bi2d_noisy = get_noisydata(bi2d_true, x, 200)
    #bi2d(bi2d_true, bi2d_noisy)


# routine for mpi parallel computations of rits
def mpi_parallel_computation():
    param_sets_sets = load_param_sets_sets()
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    psize = comm.Get_size()
    inpfile = 'rits_initial.inp.' + str(rank)
    #ns = len(param_sets_sets)
    iniidx = int(sys.argv[1])
    ns = 240
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
        bi2dt = bi2dt.reshape((ns, bi2d_true.shape[0], bi2d_true.shape[1]))
        with open('/home/kazu/desktop/240424/bi2d/bi2d_test.pkl.'+str(iniidx),
                  'wb') as f:
            pickle.dump(bi2dt, f, 4)
            pickle.dump(x, f, 4)


def gather_bi2d(timescale=600, nidx=100):
    datat = []
    for iniidx in range(nidx):
        with open('/home/kazu/desktop/240424/bi2d/bi2d_test.pkl.'
                  + str(iniidx), 'rb') as f:
            data = pickle.load(f)
            x = pickle.load(f)
        data = data[np.isnan(data).sum(axis=1).sum(axis=1) == 0]
        if iniidx == 0:
            datat = data
        else:
            datat = np.vstack((datat, data))
    datat_noisy = get_noisydata(datat, x, timescale)
    datasets = {}
    datasets['target'] = datat*timescale
    datasets['noisy'] = datat_noisy
    datasets['x'] = x
    with open('/home/kazu/desktop/240424/bi2d/bi2d_test.pkl', 'wb') as f:
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
cycles(ns=1, dim=2)
#mpi_parallel_computation()
#gather_bi2d(timescale=50, nidx=125)
#select_bi2d()
#check_data()

