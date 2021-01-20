#!/usr/bin/env python
# test of convolution NN on neutron diffraction data
# k-CV method is adopted to evaluate the robustness of the model peformance.
# each cyle in k-CV is distributed by mpi4py, so you should type mpirun -np [k-size] test2_nd_mpi.py
# 
import sys
import os
import numpy as np
sys.path.append(os.pardir)
from dataset.nd import load_nd
from conv_dp_net_nd import convdpNet
from common.optimizer import Adam
import pickle
from sklearn.model_selection import KFold
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


def run():
    num_iters = 2000
    batch_size = 100
    lossfunc = []
    acc_validation = []
    acc_train = []

    (_x_train, _t_train), (x_test, t_test) = load_nd()
    _x_train = _x_train[:,:-1]
    x_test = x_test[:,:-1]

    _ave_x_train = np.average(_x_train)
    _var_x_train = np.var(_x_train)
    _x_train = (_x_train - _ave_x_train)/np.sqrt(_var_x_train)
    x_test = (x_test - _ave_x_train)/np.sqrt(_var_x_train)

    _x_train_shape = _x_train.shape
    _x_train = _x_train.reshape(_x_train_shape[0], 1, 1, _x_train_shape[1])
    x_test_shape = x_test.shape
    x_test = x_test.reshape(x_test_shape[0], 1, 1,  x_test_shape[1])

    kf = KFold(n_splits=10, shuffle=True)
    kidx = 0
    for train_indx, validation_indx in kf.split(_x_train):
        kidx += 1
        if kidx == rank+1:
            print("validation cycle # {} starts".format(kidx))
            x_train = _x_train[train_indx]
            t_train = _t_train[train_indx]
            x_validation = _x_train[validation_indx]
            t_validation = _t_train[validation_indx]
            train_size = x_train.shape[0]

            nt = convdpNet(input_dim=(1, 1, _x_train_shape[1]),
                conv_params={'filter_num':10, 'filter_height':1, 'filter_width':5, 'pad':0, 'stride':1},
                hsize = 100, osize = 4, wstd=0.01)

            opt = Adam(lr=0.001)
            itersperepoch = max(train_size // batch_size, 1)
            for lidx in range(num_iters):
                batch_mask = np.random.choice(train_size, batch_size)
                x_batch = x_train[batch_mask]
                t_batch = t_train[batch_mask]
                #grad = nt.numerical_grad(x_batch, t_batch)
                grad = nt.analytical_grad(x_batch, t_batch)
                for key in grad.keys():
                    params = nt.params
                    opt.update(params, grad)
                #print("iteration# {} finished, lossfunc: {}".format(lidx,nt.loss(x_batch, t_batch)))
                lossfunc.append(nt.loss(x_batch, t_batch))
                if lidx % itersperepoch == 0:
                    acc_validation.append(nt.acc(x_validation, t_validation))
                    acc_train.append(nt.acc(x_train, t_train))
                    #print(nt.acc(x_test, t_test))
                    #print(nt.acc(x_train, t_train))

            print(np.average(acc_validation[-10:-1]))
    MPI.COMM_WORLD.barrier()
    all_acc_validation = comm.gather(acc_validation, root=0)
    all_acc_train  = comm.gather(acc_train, root=0)
    all_lossfunc = comm.gather(lossfunc, root=0)
    if rank == 0:
        ave_acc_validation = np.average(np.array(all_acc_validation)[:,-10:-1])
        std_acc_validation = np.var(np.array(all_acc_validation)[:,-10:-1])**0.5
        print("average acc in validation:{} +- {}".format(ave_acc_validation, std_acc_validation))
    
        pb = open('dump2_nd_cv.pkl', 'wb')
        #pickle.dump(nt, pb)
        pickle.dump(all_lossfunc, pb)
        pickle.dump(all_acc_validation, pb)
        pickle.dump(all_acc_train, pb)
        pb.close()


run()
