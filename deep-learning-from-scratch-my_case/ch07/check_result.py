#!/usr/bin/env python
#coding: utf-8
import pickle
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(6, 9))
fig.suptitle("Summary of test2_nd.py, convolution_neural-net on nd data")

def run():
    pb = open('dump2_nd_bn3_modified.pkl', 'rb')
    nt = pickle.load(pb)
    lossfunc = pickle.load(pb)
    acc_test = pickle.load(pb)
    acc_train = pickle.load(pb)
    pb.close()

    ax = fig.add_subplot(2, 1, 1)
    ax.plot(lossfunc)
    ax.set_xlabel('# of iterations')
    ax.set_ylabel('cross entropy error')
    ax = fig.add_subplot(2, 1, 2)
    ax.plot(acc_test, linewidth=2, label="acc_test")
    ax.plot(acc_train, linewidth=1, label="acc_train")
    ax.set_xlabel('# of epochs')
    ax.set_ylabel('accuracy')
    ax.set_ylim(0,1.1)
    plt.legend(loc='best')
    plt.savefig("summary_test2_nd_bn3_modified.pdf")
    plt.show()


run()
