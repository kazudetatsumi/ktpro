#!/usr/bin/env python
#coding: utf-8
import pickle
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(6, 9))
fig.suptitle("Summary of test2_nd.py, convolution_neural-net on nd data")

def run():
    pb = open('dump2_nd_bn3_modified_cv.pkl', 'rb')
    #nt = pickle.load(pb)
    lossfunc = np.array(pickle.load(pb))
    acc_test = np.array(pickle.load(pb))
    acc_train = np.array(pickle.load(pb))
    pb.close()

    ax = fig.add_subplot(2, 1, 1)
    ax.plot(np.average(lossfunc, axis=0))
    ax.set_xlabel('# of iterations')
    ax.set_ylabel('cross entropy error')
    ax = fig.add_subplot(2, 1, 2)
    #for kidx in range(0,10):
    #    ax.plot(acc_test[kidx][:])
    ave_acc_test = np.average(acc_test, axis=0)
    ave_acc_train = np.average(acc_train, axis=0)
    var_acc_test = np.var(acc_test, axis=0)
    ax.plot(ave_acc_test)
    ax.plot(ave_acc_train)
    ax.plot(ave_acc_test + var_acc_test**0.5)
    ax.plot(ave_acc_test - var_acc_test**0.5)
    #ax.errorbar(range(*ave_acc_test.shape),ave_acc_test, yerr=var_acc_test**0.5, capsize=2, fmt='o', markersize=2, ecolor='k', markeredgecolor = "k", color='none', label="acc_test")


    #ax.plot(np.average(np.array(acc_test), axis=0), linewidth=2, label="acc_test")

    #ax.plot(acc_train, linewidth=1, label="acc_train")
    ax.set_xlabel('# of epochs')
    ax.set_ylabel('accuracy')
    ax.set_ylim(0,1.1)
    plt.legend(loc='best')
    plt.savefig("summary_test2_nd.pdf")
    plt.show()


run()
