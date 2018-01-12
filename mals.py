#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

#from PIL import Image
from skimage import io

numcom = 2

def run():
   data = io.imread("Signal SI32.tif")
   datasize = data.shape
   matsize = (datasize[0],datasize[1]*datasize[2])
   x = np.reshape(data,matsize)
   print x.shape
   b = np.random.rand(numcom,matsize[0])
   print np.dot(b,x).shape
   #a = np.dot(np.dot(x,b.T),np.linalg.inv(np.dot(b,b.T)))
   #plt.plot(data1d[:,0])
   #plt.show()
run()
