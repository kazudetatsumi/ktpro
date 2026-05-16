#!/usr/bin/env python
import sys
import os
import numpy as np
sys.path.append(os.pardir)
from dataset.nd import load_nd
from conv_net_nd import convNet
from common.optimizer import Adam
import pickle


def run():
    num_iters = 10000
    batch_size = 100
    lossfunc = []
    acc_test = []
    acc_train = []

    (x_train, t_train), (x_test, t_test) = load_nd()
    x_train = x_train[:,:-1]
    x_test = x_test[:,:-1]
    ave_x_train = np.average(x_train)
    var_x_train = np.var(x_train)
    x_train = (x_train - ave_x_train)/np.sqrt(var_x_train)
    x_test = (x_test - ave_x_train)/np.sqrt(var_x_train)

    x_train_shape = x_train.shape
    x_train = x_train.reshape(x_train_shape[0], 1, 1,  x_train_shape[1])
    x_test_shape = x_test.shape
    x_test = x_test.reshape(x_test_shape[0], 1, 1,  x_test_shape[1])
    train_size = x_train_shape[0]
    nt = convNet(input_dim=(1, 1, x_train_shape[1]),
            conv_params={'filter_num':30, 'filter_height':1, 'filter_width':3, 'pad':0, 'stride':1},
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
            acc_test.append(nt.acc(x_test, t_test))
            acc_train.append(nt.acc(x_train, t_train))
            #print(nt.acc(x_test, t_test))
            #print(nt.acc(x_train, t_train))

    pb = open('dump2_nd_30_3_100.pkl', 'wb')
    pickle.dump(nt, pb)
    pickle.dump(lossfunc, pb)
    pickle.dump(acc_test, pb)
    pickle.dump(acc_train, pb)
    pb.close()


run()
