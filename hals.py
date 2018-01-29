#!/usr/bin/env python
import numpy as np
import dm4reader
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy.interpolate import griddata
import matplotlib
import sys
from libnmf import NMF

#matplotlib.font_manager._rebuild()
args = sys.argv
est=int(round((375-175)*2.5))
efi=int(round((650-225)*2.5))
numcom=2
maxcycle=1000
datanum=14
IsAligned='true'
Is3ddata='true'

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
    np_array = np_array.transpose(2,1,0) # to fit the data shape to libnmf
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
    np_array = np_array.T # to fit the data shape to libnmf
    return(np_array)


def get_1ddata(filename):
    dm4data = dm4reader.DM4File.open(filename)
    tags = dm4data.read_directory()
    image_data_tag = tags.named_subdirs['ImageList'].unnamed_subdirs[1].named_subdirs['ImageData']
    image_tag = image_data_tag.named_tags['Data']

    XDim = dm4data.read_tag_data(image_data_tag.named_subdirs['Dimensions'].unnamed_tags[0])

    np_array = np.array(dm4data.read_tag_data(image_tag))
    return(np_array)

def plotter(s):
    n=numcom
    plt.subplot2grid((4,n+3),(0,0),rowspan=1, colspan=2)
    for i in range(0,n):
       labelname="comp. #"+str(i)
       plt.plot(s[:,i],label=labelname)
    plt.tick_params(which='both',tickdir='in',top='true',right='true')
    plt.title("spectra")
    plt.legend()

def plotter2d(c,a):
    n=numcom
    for i in range(0,n):
       plt.subplot2grid((4,n+3),(0,i+2),rowspan=1, colspan=1)
       labelname="comp. #"+str(i)
       plt.imshow(c[i,:,:], cmap = cm.gray, interpolation = 'none',label="comp. #"+str(i))
       plt.title("map of comp. #"+str(i))
       plt.tick_params(labelbottom='off',bottom='false',labelleft='off',left='false')
    plt.subplot2grid((4,n+3),(0,n+2),rowspan=1, colspan=1)
    plt.imshow(a, cmap = cm.gray, interpolation = 'none',label="comp. #"+str(i))
    plt.tick_params(labelbottom='off',bottom='false',labelleft='off',left='false')
    plt.title("ADF image")
    
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
    plt.subplot2grid((4,n+3),(1,2),rowspan=3, colspan=3)
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
        n_ch = np.arange(est,efi)
        ximage = ximage[:,:,n_ch]
        xdim, ydim, Nch = ximage.shape
        X = np.reshape(ximage, (xdim*ydim, Nch))
        scale_X = np.mean(X)
        X = X / scale_X
        
        nmf = NMF(n_components=numcom, reps=3, max_itr=maxcycle)
        nmf.fit(X, num_xy=(xdim,ydim), channel_vals=n_ch, unit_name='Channel') 
        adfimage = get_2ddata(h)
        plotter(nmf.S_)
        plotter2d(nmf.C_.T.reshape([numcom,xdim,ydim]).transpose(0,2,1),adfimage)
    else:
        ximage = get_2ddata(f).T
        n_ch = np.arange(est,efi)
        ximage = ximage[:,n_ch]
        xdim, Nch = ximage.shape
        ydim = 1
        X = ximage
        scale_X = np.mean(X)
        X = X / scale_X
        nmf = NMF(n_components=numcom, reps=3, max_itr=maxcycle)
        nmf.fit(X, num_xy=(xdim,ydim), channel_vals=n_ch, unit_name='Channel') 

        adfimage = get_1ddata(h)
        plotter(nmf.S_)
        plotter1d(nmf.C_.T,adfimage)
    plotteradfsurvey(g,numcom)
    plt.show()
    #plt.savefig(outname+".eps")

run()
