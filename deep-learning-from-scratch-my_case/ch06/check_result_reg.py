#!/usr/bin/env python
# coding: utf-8
import pickle
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(6, 12))
fig.suptitle("Summary of test3.py, three-layer-neural-net on mnist data")

def run():
    pb = open('dump_reg.pkl', 'rb')
    nt = pickle.load(pb)
    lossfunc = pickle.load(pb)
    loss_test = pickle.load(pb)
    loss_train = pickle.load(pb)
    #x_train = pickle.load(pb)
    #y_train = pickle.load(pb)
    #x_test = pickle.load(pb)
    #y_test = pickle.load(pb)

    pb.close()

    ax = fig.add_subplot(4, 1, 1)
    ax.plot(lossfunc)
    ax.set_xlabel('# of iterations')
    ax.set_ylabel('mean squared error')
    ax = fig.add_subplot(4, 1, 2)
    ax.plot(loss_test, linewidth=2, label="test")
    ax.plot(loss_train, linewidth=1, label="train")
    ax.set_xlabel('# of epochs')
    ax.set_ylabel('mean squared error')
    #ax.set_ylim(0,1.1)
    ##plt.legend(loc='best')
    #ax = fig.add_subplot(4, 1, 3)
    #ax.scatter(x_train, y_train)
    #ax = fig.add_subplot(4, 1, 4)
    #ax.scatter(x_test, y_test)
    plt.savefig("summary_test_reg.pdf")
    plt.show()


run()
