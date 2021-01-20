#!/usr/bin/env python
import sys
import os
import numpy as np
sys.path.append(os.pardir)
from dataset.nd import load_nd
#from conv_bn_net_nd import convNetbn
from conv_net_ndK2 import convNetK2
from common.optimizer import Adam, SGD, Adagrad
import pickle


def run():
    np.random.seed(1671) # for reproducibility
    num_iters = 1000
    batch_size = 128
    loss_train = []
    loss_val = []
    acc_train = []
    acc_val = []

    (x_train, t_train), (x_test, t_test) = load_nd()
    xave = np.average(x_train)
    xvar = np.var(x_train)
    x_train = (x_train - xave) / np.sqrt(xvar)
    x_test = (x_test - xave) / np.sqrt(xvar)
    x_train = x_train.reshape(x_train.shape[0], 1, 1, x_train.shape[1])
    x_test = x_test.reshape(x_test.shape[0], 1, 1, x_test.shape[1])
    test_mask = np.random.choice(x_train.shape[0], x_train.shape[0]//10)
    x_validation = x_train[test_mask]
    t_validation = t_train[test_mask]
    x_train = x_train[np.invert(test_mask)]
    t_train = t_train[np.invert(test_mask)]

    nt = convNetK2(input_dim=(1, 1, x_train.shape[3]),
                   conv_params1={'FN':20, 'FH':1, 'FW':3, 'P':1, 'S':1},
                   pool_params1={'PH':1, 'PW':2, 'P':0, 'S':2},
                   conv_params2={'FN':50, 'FH':1, 'FW':5, 'P':2, 'S':1},
                   pool_params2={'PH':1, 'PW':2, 'P':0, 'S':2},
                   hsize=500, osize=4, wstd=0.01)

    opt = Adagrad(learning_rate=0.001)
    itersperepoch = max(x_train.shape[0] // batch_size, 1)
    for lidx in range(num_iters):
        batch_mask = np.random.choice(x_train.shape[0], batch_size)
        x_batch = x_train[batch_mask]
        t_batch = t_train[batch_mask]
        grad = nt.analytical_grad(x_batch, t_batch)
        for key in grad.keys():
            params = nt.params
            opt.update(params, grad)
        if lidx % itersperepoch == 0:
            loss_train.append(nt.loss(x_train, t_train))
            acc_train.append(nt.acc(x_train, t_train))
            loss_val.append(nt.loss(x_validation, t_validation))
            acc_val.append(nt.acc(x_validation, t_validation))
            print("epoch# {} ".format(lidx // itersperepoch) +
                  "loss_train: {0:.4f}, ".format(loss_train[-1]) +
                  "acc_train: {0:.4f}, ".format(acc_train[-1]) +
                  "loss_val: {0:.4f}, acc_val: {1:.4f}"
                  .format(loss_val[-1], acc_val[-1]))

    print(np.average(nt.acc(x_test, t_test)))

    pb = open('dump2_nd.pkl', 'wb')
    #pickle.dump(nt, pb)
    pickle.dump(loss_train, pb)
    pickle.dump(loss_val, pb)
    pickle.dump(acc_train, pb)
    pickle.dump(acc_val, pb)
    pb.close()


run()