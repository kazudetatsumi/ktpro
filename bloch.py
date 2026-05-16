#!/usr/bin/env python
#from mpmath import *
import math  as m
import numpy as np
import matplotlib.pyplot as plt
#from PIL import Image as im
#from mpl_toolkits.mplot3d import Axes3D
#from scipy import misc

a = 8.05
y = 0
Fg = 26.609
e0 = 200000
Vc = a**3
g = 4/a
d = a/4
w = -0.5
#beta = acot(w)
beta = 1.57 - m.atan(w)
print beta
#l = 12.26 / ((e0*(1+0.0000009788*e0))**0.5)
#theta = asin(l/(2*d))
#xai = pi*Vc*cos(theta)/(l*Fg)
#s = w/xai
#kk = 1/l
#kkz = kk*cos(2*theta)
#kx = kk*sin(2*theta)
##ky = 0
#y = 0
#k1z = kkz + (s - (s**2 + 1/xai**2)**0.5)/2
#k2z = kkz + (s + (s**2 + 1/xai**2)**0.5)/2
#x = np.linspace(0, 16.1, 100)
#z = np.linspace(0, 2500, 100)
#b1 = (2*pi*j*(kx*x+ky*y+k1z*z))
#print exp(b1[0])
