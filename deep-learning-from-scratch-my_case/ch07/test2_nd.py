#!/usr/bin/env python
import sys
import os
import numpy as np
sys.path.append(os.pardir)
from dataset.nd import load_nd
#from conv_bn_net_nd import convNetbn
from conv_net_nd import convNet
from common.optimizer import Adam
import pickle
from sklearn.model_selection import KFold


def run():
    num_iters = 1000
    batch_size = 40
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
        if kidx == 1:
            print("validation cycle # {} starts".format(kidx))
            x_train = _x_train[train_indx]
            t_train = _t_train[train_indx]
            x_validation = _x_train[validation_indx]
            t_validation = _t_train[validation_indx]
            train_size = x_train.shape[0]

            nt = convNet(input_dim=(1, 1, _x_train_shape[1]),
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
                print("iteration# {} finished, lossfunc: {}".format(lidx,nt.loss(x_batch, t_batch)))
                lossfunc.append(nt.loss(x_batch, t_batch))
                if lidx % itersperepoch == 0:
                    acc_validation.append(nt.acc(x_validation, t_validation))
                    acc_train.append(nt.acc(x_train, t_train))
                    #print(nt.acc(x_test, t_test))
                    #print(nt.acc(x_train, t_train))

            print(np.average(acc_validation[-50:-1]))

        pb = open('dump2_nd.pkl', 'wb')
        pickle.dump(nt, pb)
        pickle.dump(lossfunc, pb)
        pickle.dump(acc_validation, pb)
        pickle.dump(acc_train, pb)
        pb.close()


run()
