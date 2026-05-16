#!/usr/bin/env python
import sys
sys.path.append("/home/kazu/ktpro")
import count_unmasked_elements as cue

def run():
    hdffile = "/home/kazu/desktop/200204/fine/hourbyhour/" + \
              "ortho_opt_without_mask/10h/eliminated_data.hdf5"
    pro = cue.Count_Unmasked_Elements(hdffile)
    print("num of unmasked elements:", pro.get_size())
    print("frac of unmasked elements:", pro.get_frac())

    pro.file = "/home/kazu/desktop/200522/Ei42/veryfineq/14m/eliminated_data.hdf5"
    pro.get_arrays()
    print("num of unmasked elements:", pro.get_size())
    print("frac of unmasked elements:",
          pro.get_frac())

    pro.file = "/home/kazu/desktop/200522/Ei24/fineq/26m/eliminated_data.hdf5"
    pro.get_arrays()
    print("num of unmasked elements:", pro.get_size())
    print("frac of unmasked elements:", pro.get_frac())



run()
