#!/usr/bin/env python
import numpy as np
import dm4reader
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.interpolate import griddata
import matplotlib
import sys

#matplotlib.font_manager._rebuild()
args = sys.argv
est=int(round((375-175)*2.5))
efi=int(round((650-225)*2.5))
#est=int(round((375-175)))
#efi=int(round((650-225)))
numcom=2
maxcycle=100
datanum=5
#datanum=args[1]
IsAligned='true'
Is3ddata='true'
#IsAligned=args[2]

if IsAligned == 'true':
    f="/home/kazu/LaNi5_dm4/SI data ("+str(datanum)+")Signal SI_spike_bkg_rmvd (aligned).dm4"
    outname="result_no"+str(datanum)+"_aligned_comp"+str(numcom)+"_cycle"+str(maxcycle)
else:
    f="/home/kazu/LaNi5_dm4/SI data ("+str(datanum)+")Signal SI_spike_bkg_rmvd.dm4"
    outname="result_no"+str(datanum)+"_comp"+str(numcom)+"_cycle"+str(maxcycle)
#g="/home/kazu/LaNi5_dm4/ADF Image (SI Survey)"+str(datanum)+".tif"
g="/home/kazu/LaNi5_dm4/ADF Image"+str(datanum)+" (SI Survey).tif"
h="/home/kazu/LaNi5_dm4/ADF Image"+str(datanum)+".dm4"



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

def get_2ddata(filename):
    dm4data = dm4reader.DM4File.open(filename)
    tags = dm4data.read_directory()
    image_data_tag = tags.named_subdirs['ImageList'].unnamed_subdirs[1].named_subdirs['ImageData']
    image_tag = image_data_tag.named_tags['Data']

    XDim = dm4data.read_tag_data(image_data_tag.named_subdirs['Dimensions'].unnamed_tags[0])
    YDim = dm4data.read_tag_data(image_data_tag.named_subdirs['Dimensions'].unnamed_tags[1])

    np_array = np.array(dm4data.read_tag_data(image_tag))
    np_array = np.reshape(np_array, (YDim, XDim))
    return(np_array)


def get_1ddata(filename):
    dm4data = dm4reader.DM4File.open(filename)
    tags = dm4data.read_directory()
    image_data_tag = tags.named_subdirs['ImageList'].unnamed_subdirs[1].named_subdirs['ImageData']
    image_tag = image_data_tag.named_tags['Data']

    XDim = dm4data.read_tag_data(image_data_tag.named_subdirs['Dimensions'].unnamed_tags[0])

    np_array = np.array(dm4data.read_tag_data(image_tag))
    return(np_array)

def mals3d(data,n,mc):
    srs=[]
    dim = data.shape
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
        if nn > maxcycle*0.2:
            print srs[nn-1] - sr
        nn = nn+1
    print sr
    b2d = np.reshape(b,[n,dim[1],dim[2]])
    print b2d.shape
    return(a,b2d,srs)

def mals2d(data,n,mc):
    srs=[]
    ddim = data.shape
    x = data
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
        if nn > maxcycle*0.2:
            print srs[nn-1] - sr
        nn = nn+1
    print sr
    return(a,b,srs)

def plotter(s):
    n=numcom
    #subplot2grid(shape, loc, rowspan=1, colspan=1, fig=None, **kwargs
    #ax1 = plt.subplot2grid((6,1), (0,0), rowspan=1, colspan=1)
    #ax2 = plt.subplot2grid((6,1), (1,0), rowspan=4, colspan=1)
    #ax3 = plt.subplot2grid((6,1), (5,0), rowspan=1, colspan=1)
    #plt.subplot(1,n+2,1)
    plt.subplot2grid((4,n+3),(0,0),rowspan=1, colspan=2)
    for i in range(0,n):
       labelname="comp. #"+str(i)
       plt.plot(s[:,i],label=labelname)
    #plt.plot(s[:,1],label="comp. #2")
    plt.tick_params(which='both',tickdir='in',top='true',right='true')
    plt.title("spectra")
    plt.legend()

def plotter2d(c,a):
    n=numcom
    for i in range(0,n):
    #   plt.subplot(1,n+2,i+2)
    #plt.subplot2grid((1,4),(0,2),rowspan=1, colspan=1)
       plt.subplot2grid((4,n+3),(0,i+2),rowspan=1, colspan=1)
       labelname="comp. #"+str(i)
       plt.imshow(c[i,:,:], cmap = cm.gray, interpolation = 'none',label="comp. #"+str(i))
       plt.title("map of comp. #"+str(i))
       plt.tick_params(labelbottom='off',bottom='false',labelleft='off',left='false')
    plt.subplot2grid((4,n+3),(0,n+2),rowspan=1, colspan=1)
    plt.imshow(a, cmap = cm.gray, interpolation = 'none',label="comp. #"+str(i))
    plt.tick_params(labelbottom='off',bottom='false',labelleft='off',left='false')
    plt.title("ADF image")
    #plt.subplot2grid((1,4),(0,3),rowspan=1, colspan=1)
    #plt.subplot(1,3,3)
    #plt.imshow(c[1,:,:], cmap = cm.gray, interpolation = 'none',label="comp. #2")
    #plt.title("map of comp. #2")
    
def plotter1d(c,a):
    n=numcom
    for i in range(0,n):
       plt.subplot2grid((4,n+3),(0,i+2),rowspan=1, colspan=1)
       labelname="comp. #"+str(i)
       plt.plot(c[i,:], label="comp. #"+str(i))
       plt.title("map of comp. #"+str(i))
       plt.tick_params(labelbottom='off',bottom='false',labelleft='off',left='false')
    plt.subplot2grid((4,n+3),(0,n+2),rowspan=1, colspan=1)
    plt.plot(a, label="adf")
    plt.tick_params(labelbottom='off',bottom='false',labelleft='off',left='false')
    plt.title("ADF image")
    
def plotteradfsurvey(tifname,n):
    n=numcom
    adfimage=plt.imread(tifname)
    #plt.subplot2grid((4,n+3),(0,n+2),rowspan=1, colspan=1)
    #plt.tick_params(labelbottom='off',bottom='false',labelleft='off',left='false')
    #plt.title("ADF survey image")
    #plt.imshow(adfimage)
    plt.subplot2grid((4,n+3),(1,2),rowspan=3, colspan=3)
    #plt.subplot(1,n+2,n+2)
    plt.imshow(adfimage)
    if Is3ddata == 'true':
       plt.subplots_adjust(top=0.95,hspace=0.01,wspace=0.01)
    else:
       plt.subplots_adjust(top=0.93,hspace=0.10,wspace=0.01)
    plt.tick_params(labelbottom='off',bottom='false',labelleft='off',left='false')
    plt.title("ADF survey image")
    


def run():
    plt.figure(figsize=(2.5*(numcom+3),2.5*4))
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.suptitle(outname)
    if Is3ddata == 'true':
        ximage = get_3ddata(f)
        ximage2 = ximage[est:efi,:,:]
        adat,bdat,sumress = mals3d(ximage2,numcom,maxcycle)
        adfimage = get_2ddata(h)
        plotter(adat)
        plotter2d(bdat,adfimage)
    else:
        ximage = get_2ddata(f).T
        ximage2 = ximage[est:efi,:]
        adat,bdat,sumress = mals2d(ximage2,numcom,maxcycle)
        adfimage = get_1ddata(h)
        plotter(adat)
        plotter1d(bdat,adfimage)
    plotteradfsurvey(g,numcom)
    #plt.show()
    plt.savefig(outname+".eps")

run()
