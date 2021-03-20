#!/usr/bin/env python
import numpy as np
import math
import matplotlib.pyplot as plt


class GaussianProcess:
    
    def __init__(self, num_inputs, min_x=-5, max_x=5, ms=".", ls="", c="k"):
        self.num_inputs = num_inputs
        self.min_x = min_x
        self.max_x = max_x
        self.ms = ms
        self.ls = ls
        self.c = c
        self.generate_input()
        self.create_fig()

    def create_fig(self, title="Gaussian Process exercise"):
        self.fig = plt.figure(figsize=(6, 4.5))
        self.fig.suptitle(title)

    def generate_input(self):
        self.x = np.linspace(self.min_x, self.max_x, self.num_inputs)

    def generate_L(self):
        K = np.zeros((self.num_inputs, self.num_inputs))
        for p, xp in enumerate(self.x):
            for q, xq in enumerate(self.x):
                K[p, q] = math.exp(-0.5*abs(xp - xq)**2)
        K += 0.00000001*np.identity(self.num_inputs)
        self.L = np.linalg.cholesky(K)

    def output(self):
        u = np.random.normal(size=self.num_inputs)
        f = np.dot(self.L, u)

        plt.plot(self.x, f, c=self.c, marker=self.ms, linestyle=self.ls)


def run():
    num_inputs = 80
    prj = GaussianProcess(num_inputs)
    prj.generate_L()
    prj.output()
    prj.ms = ""
    prj.ls = "-"
    prj.c = "r"
    prj.output()
    prj.c = "g"
    prj.output()
    plt.xlabel("input, x")
    plt.ylabel("output, f(x)")


run()
plt.show()



