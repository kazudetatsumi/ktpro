#!/usr/bin/env python
import numpy as np
import dm4reader
import matplotlib.pyplot as plt
import math

f="tmp.dm4"
yst=1
xst=1
est=int(round((375-175)*2.5))
efi=int(round((650-225)*2.5))
numcom=2
maxcycle=2000


def get_3ddata(filename):
    dm4data = dm4reader.DM4File.open("tmp.dm4")
    tags = dm4data.read_directory()
    image_data_tag = tags.named_subdirs['ImageList'].unnamed_subdirs[1].named_subdirs['ImageData']
    image_tag = image_data_tag.named_tags['Data']

    XDim = dm4data.read_tag_data(image_data_tag.named_subdirs['Dimensions'].unnamed_tags[0])
    YDim = dm4data.read_tag_data(image_data_tag.named_subdirs['Dimensions'].unnamed_tags[1])
    ZDim = dm4data.read_tag_data(image_data_tag.named_subdirs['Dimensions'].unnamed_tags[2])

    np_array = np.array(dm4data.read_tag_data(image_tag))
    np_array = np.reshape(np_array, (ZDim, YDim, XDim))
    return(np_array)

def mals(data,n,mc):
    srs=[]
    dim = data.shape
    xdim = dim[2]
    ydim = dim[1]
    ddim = [dim[0],dim[2]*dim[1]]
    X = np.reshape(data,ddim)
    b = np.random.rand(n,ddim[1])
    a = np.dot( np.dot( x, b.T ), np.linalg.inv( np.dot( b, b.T ) ) )
    nn = 0
    while nn < mc:
        for i in range(0,n):
            a[:,i] = a[:,i] / np.linalg.norm(a[:,i])
        a = np.where( a > 0, a, 0)
        wb = np.eye(numcom)*math.exp(-(float(nn/100)))
        b = np.dot( np.linalg.inv( np.dot( a.T, a) + wb ),  np.dot( a.T, x ) + np.dot( wb, b )  )
        b = np.where( b > 0, b, 0)
        wa = np.eye(numcom)*math.exp(-(float(nn/100)))
        a = np.dot(  np.dot( x, b.T ) + np.dot( a, wa ), np.linalg.inv( np.dot( b, b.T ) + wa )  )
        r = x - np.dot(a,b)
        sr = np.linalg.norm(r)
        srs.append(sr)
        if nn > 1:
            print srs[nn-1] - sr
        nn = nn+1
    print sr
    return(a,b,srs)

def plotter(s):
    plt.plot(s[:,0])
    plt.plot(s[:,1])
    plt.show()


def run():
    ximage = get_3ddata(f)
    ximage2 = ximage[est:efi,:,:]
    adat,bdat,sumress = mals(ximage2,numcom,maxcycle)
    plotter(adat)

run()
