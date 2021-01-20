#!/usr/bin/env python
import sys
import os
sys.path.append(os.pardir)
from dataset.mnist import load_mnist
import numpy as np
import matplotlib.pyplot as plt
import pickle


def get_data():
    (x_train, t_train), (x_test, t_test) = load_mnist(
                                                     flatten=True,
                                                     normalize=True,
                                                     one_hot_label=False
                                                     )
    #print(x_train.shape)
    #print(t_train.shape)
    #print(x_test.shape)
    #print(t_test.shape)
    #print(t_train[0:100])
    #x_train = np.reshape(x_train, (60000, 28, 28))
    #plt.imshow(x_train[2, :, :])
    #print(t_train[2])
    return x_test, t_test


def get_network():
    with open("sample_weight.pkl", 'rb') as f:
        network = pickle.load(f)
    return network


def sigmoid(x):
    return 1. / (1. + np.exp(-x))


def softmax(x):
    return np.exp(x) / np.sum(np.exp(x))


def cross_entropy_err(y, t):
    if y.ndim == 1:
        t = t.reshape(1, t.size)
        y = y.reshape(1, y.size)
    batch_size = y.shape[0]
    return -np.sum(np.log(y[np.arange(batch_size), t] + 1e-7)) / batch_size


def mse(y, t):
    return 0.5*np.sum((y-t)**2)


def test(network, x):
    W1, W2, W3 = network['W1'], network['W2'], network['W3']
    b1, b2, b3 = network['b1'], network['b2'], network['b3']
    #print("size of Weight 1 matrix: {}".format(W1.shape))
    #print("size of Weight 2 matrix: {}".format(W2.shape))
    #print("size of Weight 3 matrix: {}".format(W3.shape))
    #print("size of bias 1 matrix: {}".format(b1.shape))
    #print("size of bias 2 matrix: {}".format(b2.shape))
    #print("size of bias 3 matrix: {}".format(b3.shape))

    a1 = np.dot(x, W1) + b1
    z1 = sigmoid(a1)
    a2 = np.dot(z1, W2) + b2
    z2 = sigmoid(a2)
    a3 = np.dot(z2, W3) + b3
    y = softmax(a3)
    return y




def run():
    x, t = get_data()
    print("size of data matrix: {}".format(x.shape))
    print("size of label matrix: {}". format(t.shape))
    nt = get_network()
    y = test(nt, x)
    print("size of infered label matrix: {}".format(y.shape))
    truecount = 0
    et = np.argmax(y, axis=1)
    truecount = np.sum(et == t)
    print("Fraction of samples having the infered label as the true label: {}".format(truecount/len(t)))


def run_batch():
    batch_size = 100
    x, t = get_data()
    nt = get_network()
    truecount = 0
    for sidx in range(0, x.shape[0], batch_size):
        x_batch = x[sidx:sidx+batch_size, :]
        t_batch = t[sidx:sidx+batch_size]
        y_batch = test(nt, x_batch)
        et = np.argmax(y_batch, axis=1)
        print("Cross entropy error: {}".format(cross_entropy_err(y_batch, t_batch)))
        truecount += np.sum(et == t[sidx:sidx+batch_size])
    print("Fraction of samples having the infered label as the true label: {}".format(truecount/len(t)))






run_batch()
