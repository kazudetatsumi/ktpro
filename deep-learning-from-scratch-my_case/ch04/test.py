import sys, os
sys.path.append(os.pardir)
import numpy as np
from common.functions import softmax, cross_entropy_error
from common.gradient import numerical_gradient


class Net:
    def __init__(self):
        self.W  = np.random.randn(2,3)


    def get_a1(self, x):
        return np.dot(x, self.W)


    def get_cee(self, x, t):
        a1 = self.get_a1(x)
        y = softmax(a1)
        cee = cross_entropy_error(y, t)
        return(cee)


def run():
    mynet = Net()
    x_data = np.array([0.6, 0.9])
    true_label = np.array([0, 0, 1])
    print("present weight matrix:\n {}".format(mynet.W))
    #def f(dummy):                                              # here f has a dummy argument, but f depends on mynet.W through calss Net.
    #    return mynet.get_cee(x_data, true_label)               # f does not depend on the dummy argument but mynet.W.
    #                                                           # Due to the process in numerical_gradient,  f must have one argument.
    f = lambda dummy: mynet.get_cee(x_data, true_label)         # f can be written by using a lambda equation.

    print("gradient of present wieht matrix:\n {}".format(numerical_gradient(f, mynet.W)))


run()
