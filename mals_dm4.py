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
sumress=[]


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

def run():
    ximage = get_3ddata(f)
    ximage2 = ximage[est:efi,:,:]
    dim = ximage2.shape
    xdim = dim[2]
    ydim = dim[1]
    ddim = [dim[0],dim[1]*dim[2]]
    xdat = np.reshape(ximage2,ddim)
    bdat = np.random.rand(numcom,ddim[1])
    print np.dot( np.dot( xdat, bdat.T ), np.dot( bdat, bdat.T ) ).shape
    adat = np.dot( np.dot( xdat, bdat.T ), np.linalg.inv( np.dot( bdat, bdat.T ) ) )

    
    nn=0
    while nn < maxcycle:
        for i in range(0,numcom):
            adat[:,i] = adat[:,i] / np.linalg.norm(adat[:,i])
        adat = np.where( adat > 0, adat, 0)
        wb = np.eye(numcom)*math.exp(-(float(nn/100)))
        bdat = np.dot( np.linalg.inv( np.dot( adat.T, adat) + wb ),  np.dot( adat.T, xdat ) + np.dot( wb, bdat )  )
        bdat = np.where( bdat > 0, bdat, 0)
        wa = np.eye(numcom)*math.exp(-(float(nn/100)))
        adat = np.dot(  np.dot( xdat, bdat.T ) + np.dot( adat, wa ), np.linalg.inv( np.dot( bdat, bdat.T ) + wa )  )
        res = xdat - np.dot(adat,bdat)
        sumres = np.linalg.norm(res)
        sumress.append(sumres)
        if nn > 1:
            print sumress[nn-1] - sumres
        nn = nn+1
    print sumres
    plt.plot(adat[:,0])
    plt.plot(adat[:,1])
    plt.show()

run()
