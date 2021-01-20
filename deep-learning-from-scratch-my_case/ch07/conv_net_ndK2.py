#!/usr/bin/env python
import sys
import os
sys.path.append(os.pardir)
import numpy as np
from collections import OrderedDict
from common.layers import *


class convNetK2:
    def __init__(self, input_dim=(1, 1, 2299),
                 conv_params1={'FN':30, 'FH':1, 'FW':5, 'P':0, 'S':1},
                 pool_params1={'PH':1, 'PW':5, 'P':0, 'S':1},
                 conv_params2={'FN':30, 'FH':1, 'FW':5, 'P':0, 'S':1},
                 pool_params2={'PH':1, 'PW':5, 'P':0, 'S':1},
                 hsize=100, osize=4, wstd=0.01):

        # caluculating the sizes of matrices..
        FW1 = conv_params1['FW']
        P1 = conv_params1['P']
        S1 = conv_params1['S']
        W1 = input_dim[2]
        OW1 = (W1 - FW1 + 2*P1)/S1 + 1
        PW1 = pool_params1['PW']
        PP1 = pool_params1['P']
        POW1 = int((OW1 + 2*PP1 - PW1)/PW1 + 1)

        FW2 = conv_params2['FW']
        P2 = conv_params2['P']
        S2 = conv_params2['S']
        OW2 = (POW1 - FW2 + 2*P2)/S2 + 1
        PW2 = pool_params2['PW']
        PP2 = pool_params2['P']
        POW2 = int((OW2 + 2*PP2 - PW2)/PW2 + 1)

        FN1 = conv_params1['FN']
        FH1 = conv_params1['FH']
        C1 = input_dim[0]

        FN2 = conv_params2['FN']
        FH2 = conv_params2['FH']

        pool_output_size = POW2*FN2
        print("pool_output_size",pool_output_size)

        # parameters
        self.params = {}
        #wstd = np.sqrt(2.0/(1*FW1)).astype(np.float32)
        wstd = np.sqrt(2.0/(1*FW1+FN1*FW1/2)).astype(np.float32)
        #self.params['W1'] = (wstd*np.random.randn(FN1, C1, FH1, FW1))\
        #    .astype(np.float32)
        self.params['W1'] = (np.random.uniform(
                                               -wstd, wstd,
                                               (FN1, C1, FH1, FW1)))\
            .astype(np.float32)
        self.params['b1'] = np.zeros(FN1, dtype='float32')
        #wstd = np.sqrt(2.0/(FN1*FW2)).astype(np.float32)
        wstd = np.sqrt(2.0/(FN1*FW2+FN2*FW2/2)).astype(np.float32)
        #self.params['W2'] = (wstd*np.random.randn(FN2, FN1, FH2, FW2))\
        #    .astype(np.float32)
        self.params['W2'] = (np.random.uniform(
                                               -wstd, wstd,
                                               (FN2, FN1, FH2, FW2)))\
            .astype(np.float32)
        self.params['b2'] = np.zeros(FN2, dtype='float32')
        #wstd = np.sqrt(2.0/pool_output_size).astype(np.float32)
        wstd = np.sqrt(2.0/(pool_output_size+hsize)).astype(np.float32)
        #self.params['W3'] = (wstd*np.random.randn(pool_output_size, hsize))\
        #    .astype(np.float32)
        self.params['W3'] = (np.random.uniform(
                                               -wstd, wstd,
                                               (pool_output_size, hsize)))\
            .astype(np.float32)
        self.params['b3'] = np.zeros(hsize, dtype='float32')
        #wstd = np.sqrt(2.0/hsize).astype(np.float32)
        wstd = np.sqrt(2.0/(hsize+osize)).astype(np.float32)
        self.params['W4'] = (np.random.uniform(
                                               -wstd, wstd,
                                               (hsize, osize)))\
            .astype(np.float32)
        #self.params['W4'] = (wstd*np.random.randn(hsize, osize))\
        #    .astype(np.float32)
        self.params['b4'] = np.zeros(osize, dtype='float32')

        # layers
        self.layers = OrderedDict()
        self.layers['Conv1'] = ConvolutionK(self.params['W1'],
                                            self.params['b1'],
                                            1, conv_params1['S'],
                                            0, conv_params1['P'])
        self.layers['Relu1'] = Relu()
        self.layers['Pool1'] = PoolingK(pool_params1['PH'],
                                        pool_params1['PW'],
                                        1, pool_params1['S'],
                                        0, pool_params1['P'])
        self.layers['Conv2'] = ConvolutionK(self.params['W2'],
                                            self.params['b2'],
                                            1, conv_params2['S'],
                                            0, conv_params2['P'])
        self.layers['Relu2'] = Relu()
        self.layers['Pool2'] = PoolingK(pool_params2['PH'],
                                        pool_params2['PW'],
                                        1, pool_params2['S'],
                                        0, pool_params2['P'])
        self.layers['Affine1'] = Affine(self.params['W3'], self.params['b3'])
        self.layers['Relu3'] = Relu()
        self.layers['Affine2'] = Affine(self.params['W4'], self.params['b4'])

        self.lastLayer = SoftmaxWithLoss()

    def predict(self, x):
        for layer in self.layers.values():
            #print("x", x.shape)
            #print(layer)
            x = layer.forward(x)
        return x

    def loss(self, x, t):
        y = self.predict(x)
        return self.lastLayer.forward(y, t)

    def acc(self, x, t):
        y = self.predict(x)
        if t.ndim != 1:
            t = np.argmax(t, axis=1)
        acc = np.sum(np.argmax(y, axis=1) == t) / float(t.shape[0])
        return acc

    def analytical_grad(self, x, t):
        # forward   Is this necessary? -> Yes. It conveys t.
        self.loss(x, t)

        # backward
        dout = 1
        dout = self.lastLayer.backward(dout)
        layers = list(self.layers.values())
        layers.reverse()
        for layer in layers:
            dout = layer.backward(dout)

        grads = {}
        grads['W1'] = self.layers['Conv1'].dW
        grads['b1'] = self.layers['Conv1'].db
        grads['W2'] = self.layers['Conv2'].dW
        grads['b2'] = self.layers['Conv2'].db
        grads['W3'] = self.layers['Affine1'].dW
        grads['b3'] = self.layers['Affine1'].db
        grads['W4'] = self.layers['Affine2'].dW
        grads['b4'] = self.layers['Affine2'].db

        return grads

