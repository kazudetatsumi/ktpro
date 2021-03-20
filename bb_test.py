#!/usr/bin/env python
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
import gp_test as gp
import matplotlib.pyplot as plt


class BrownianBridgeProcess(gp.GaussianProcess):
    def __init__(self, num_inputs=80, min_x=0.0, max_x=1.00):
        self.num_inputs = num_inputs
        self.min_x = min_x
        self.max_x = max_x
        super(BrownianBridgeProcess, self).generate_input()

    def calc_K(self, x1, x2):
        K = np.zeros((x1.shape[0], x2.shape[0]))
        for p, xp in enumerate(x1):
            for q, xq in enumerate(x2):
                K[p, q] = min(xp, xq) - xp*xq
        return(K)

    def proc(self):
        super(BrownianBridgeProcess, self).calc_L(self.x)


def samplerun():
    prj = BrownianBridgeProcess()
    prj.proc()
    for i in range(3):
        prj.draw_random_func()
    plt.xlabel("input, x")
    plt.ylabel("output, f(x)")


samplerun()
plt.show()
