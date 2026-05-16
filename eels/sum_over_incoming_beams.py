#!/usr/bin/env python
import numpy as np
import h5py


nb = 48
proplist = [1,2,3,4]
edgelist = ['e1','e2']
headname = 'co_l23'
g = './co_l23_ddscs.hdf5'

def parse_data(filename):
    data = []
    with open(filename) as f:
        for line in f:
            if line.strip() == "":
                continue
            data.append([float(x) for x in line.split()])
        return np.array(data)

def sumdata(nprop,edge):
    tmp_filename= headname + '_' +  str(1) + '_' + str(nprop) + '_' + edge
    tmp_data=parse_data(tmp_filename)
    all_data=np.zeros_like(tmp_data)
    for i in range(1,nb+1):
        filename = headname + '_' +  str(i) + '_' + str(nprop) + '_' + edge
        print filename
        each_data=parse_data(filename) 
        all_data = all_data + each_data
    return (all_data)

def contract(propl,edgel):
    tmpsize = np.shape(sumdata(propl[0],edgel[0]))
    final_data=np.zeros((len(propl),len(edgel),tmpsize[0],tmpsize[1]))
    l=0
    for j in propl:
        l += 1
        m = 0
        for k in edgel:
            m += 1
            print l,m
            final_data[l-1,m-1,:,:] = sumdata(j,k)
    return(final_data,tmpsize)

def run():
    final2_data,tmpsize2=contract(proplist,edgelist)
    print final2_data
    outfh = h5py.File(g, 'w')
    outfh.create_dataset('ddscs', data = final2_data)
    outfh.create_dataset('n_incoming_beams', data = nb)
    outfh.create_dataset('n_thicknesses', data = tmpsize2[1])
    outfh.create_dataset('n_qpoints', data = tmpsize2[0])
    outfh.create_dataset('list_property', data = proplist)
    outfh.create_dataset('list_edge', data = edgelist)
    outfh.flush()
    outfh.close()
run()
