#!/usr/bin/env python
import numpy as np
import dm4reader
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.interpolate import griddata
import matplotlib
matplotlib.font_manager._rebuild()

est=int(round((375-175)*2.5))
efi=int(round((650-225)*2.5))
numcom=2
maxcycle=1000
datanum=17
IsAligned='true'

if IsAligned == 'true':
    f="/home/kazu/LaNi5_dm4/SI data ("+str(datanum)+")Signal SI_spike_bkg_rmvd (aligned).dm4"
    outname="result_no"+str(datanum)+"_aligned_comp"+str(numcom)+"_cycle"+str(maxcycle)
else:
    f="/home/kazu/LaNi5_dm4/SI data ("+str(datanum)+")Signal SI_spike_bkg_rmvd.dm4"
    outname="result_no"+str(datanum)+"_comp"+str(numcom)+"_cycle"+str(maxcycle)
g="/home/kazu/LaNi5_dm4/ADF Image (SI Survey)"+str(datanum)+".tif"



def get_3ddata(filename):
    dm4data = dm4reader.DM4File.open(filename)
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
    x = np.reshape(data,ddim)
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
        #print sr
        srs.append(sr)
        #if nn > 1:
        #    print srs[nn-1] - sr
        nn = nn+1
    print sr
    b2d = np.reshape(b,[n,dim[1],dim[2]])
    print b2d.shape
    return(a,b2d,srs)

def plotter(s,n):
    #subplot2grid(shape, loc, rowspan=1, colspan=1, fig=None, **kwargs
    #ax1 = plt.subplot2grid((6,1), (0,0), rowspan=1, colspan=1)
    #ax2 = plt.subplot2grid((6,1), (1,0), rowspan=4, colspan=1)
    #ax3 = plt.subplot2grid((6,1), (5,0), rowspan=1, colspan=1)
    #plt.subplot(1,n+2,1)
    plt.subplot2grid((1,n+3),(0,0),rowspan=1, colspan=2)
    for i in range(0,n):
       labelname="comp. #"+str(i)
       plt.plot(s[:,i],label=labelname)
    #plt.plot(s[:,1],label="comp. #2")
    plt.tick_params(which='both',tickdir='in',top='true',right='true')
    plt.title("spectra")
    plt.legend()

def plotter2d(c,n):
    for i in range(0,n):
    #   plt.subplot(1,n+2,i+2)
    #plt.subplot2grid((1,4),(0,2),rowspan=1, colspan=1)
       plt.subplot2grid((1,n+3),(0,i+2),rowspan=1, colspan=1)
       labelname="comp. #"+str(i)
       plt.imshow(c[i,:,:], cmap = cm.gray, interpolation = 'none',label="comp. #"+str(i))
       plt.title("map of comp. #"+str(i))
       plt.tick_params(labelbottom='off',bottom='false',labelleft='off',left='false')
    #plt.subplot2grid((1,4),(0,3),rowspan=1, colspan=1)
    #plt.subplot(1,3,3)
    #plt.imshow(c[1,:,:], cmap = cm.gray, interpolation = 'none',label="comp. #2")
    #plt.title("map of comp. #2")
    
def plotteradf(tifname,n):
    adfimage=plt.imread(tifname)
    plt.subplot2grid((1,n+3),(0,n+2),rowspan=1, colspan=1)
    #plt.subplot(1,n+2,n+2)
    plt.imshow(adfimage)
    #plt.subplots_adjust(bottom=0.1,top=0.82,wspace=0.1)
    plt.subplots_adjust(top=0.8,wspace=0.1)
    plt.tick_params(labelbottom='off',bottom='false',labelleft='off',left='false')
    plt.title("ADF survey image")
    


def run():
    ximage = get_3ddata(f)
    ximage2 = ximage[est:efi,:,:]
    adat,bdat,sumress = mals(ximage2,numcom,maxcycle)
    plt.figure(figsize=(2.5*(numcom+3),2.5))
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.suptitle(outname)
    plotter(adat,numcom)
    plotter2d(bdat,numcom)
    plotteradf(g,numcom)
    #plt.show()
    plt.savefig(outname+".eps")

run()
