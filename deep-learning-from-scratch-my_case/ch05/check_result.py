#!/usr/bin/env python
# coding: utf-8
import pickle
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(6, 9))
fig.suptitle("Summary of test3.py, three-layer-neural-net on mnist data")

def run():
    pb = open('dump4.pkl', 'rb')
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
    ax.set_ylim(0,1)
    plt.legend()
    plt.savefig("summary_test4.pdf")
    plt.show()


run()
