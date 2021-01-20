#!/usr/bin/env python
import sys
import os
import numpy as np
sys.path.append(os.pardir)
from two_layer_net import twoLNet
from common.optimizer import Adagrad
from dataset.nd import load_nd
import pickle


def run():
    num_iters = 20000
    batch_size = 200
    lossfunc = []
    acc_test = []
    acc_train = []

    (x_train, t_train), (x_test, t_test) = load_nd()
    print(x_train.shape)
    print(t_train.shape)
    #x_train = x_train / np.max(x_train)
    #x_test = x_test / np.max(x_test)

    ave_x_train = np.average(x_train)
    var_x_train = np.var(x_train)
    x_train = (x_train - ave_x_train)/np.sqrt(var_x_train)
    x_test = (x_test - ave_x_train)/np.sqrt(var_x_train)
    train_size = t_train.shape[0]
    nt = twoLNet(dsize= x_train.shape[1], hsize=50, osize=t_test.shape[1], wstd=0.01)
    opt = Adagrad()
    itersperepoch = train_size // batch_size
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
    #predicted = nt.predict(x_test)
    #print(np.sum(np.argmax(predicted,axis=1) == np.argmax(t_test,axis=1)) / float(t_test.shape[0]))
    #print(acc_test[-1])
    #print(acc_train[-1])

    pb = open('dump2_nd.pkl', 'wb')
    pickle.dump(nt, pb)
    pickle.dump(lossfunc, pb)
    pickle.dump(acc_test, pb)
    pickle.dump(acc_train, pb)
    pb.close()


run()
