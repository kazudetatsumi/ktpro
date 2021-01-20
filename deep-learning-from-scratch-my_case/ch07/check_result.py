#!/usr/bin/env python
#coding: utf-8
import pickle
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(6, 9))
fig.suptitle("Summary of test2_nd.py, convolution_neural-net on nd data")


def run():
    pb = open('dump2_nd.pkl', 'rb')
    loss_train = pickle.load(pb)
    loss_val = pickle.load(pb)
    acc_train = pickle.load(pb)
    acc_val = pickle.load(pb)
    pb.close()

    ax = fig.add_subplot(2, 1, 1)
    ax.plot(loss_val, linewidth=2, label="acc_val")
    ax.plot(loss_train, linewidth=1, label="acc_train")
    ax.set_xlabel('# of iterations')
    ax.set_ylabel('cross entropy error')
    ax = fig.add_subplot(2, 1, 2)
    ax.plot(acc_val, linewidth=2, label="acc_val")
    ax.plot(acc_train, linewidth=1, label="acc_train")
    ax.set_xlabel('# of epochs')
    ax.set_ylabel('accuracy')
    ax.set_ylim(0,1.1)
    plt.legend(loc='best')
    plt.savefig("summary_test2_nd.pdf")
    plt.show()


run()
