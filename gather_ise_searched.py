#!/usr/bin/env python
# gather lowest ise in try directries.
# Kazuyoshi TATSUMI 2020.
# gather_ise_fine method is added to tune the selection of bin-widths at ISE
# minimaISE.
# Kazuyoshi TATSUMI 2021.
import sys
sys.path.append("/home/kazu/ktpro")
import check_ise_class as cic
import numpy as np
import h5py
import subprocess


class gather_ise_search(cic.SCheckiselog):
    def __init__(self, times, head, shift=None):
        self.times = times
        self.head = head
        self.shift = shift

    def gather_ise(self):
        line = subprocess.getoutput("tail --line=1 " + self.log)
        #print(line, self.log)
        if "(" in line:
            values = line.split("(")
            values2 = values[1].split(")")
        else:
            for line in open(self.log):
                if "(" in line:
                    values = line.split("(")
                    values2 = values[1].split(")")
        self.minbw = np.array(values2[0].split(","), dtype='int32')
        self.minise = float(values2[1].split(" ")[5])

    def gather_ise_fine(self):
        super(gather_ise_search, self).__init__(self.log, self.reflog)
        super(gather_ise_search, self).getiselogdata(self._shift)
        self.minbw = np.squeeze(self.bw[self.ise == np.min(self.ise)])
        self.minise = np.min(self.ise)

    def read_h5py(self, outfile, term):
        f = h5py.File(outfile, 'r')
        return np.array(f[term])

    def get_ise_from_one_dir(self):
        self.n = np.sum(self.read_h5py(self.workdir+"/eliminated_data.hdf5",
                                       "data4")*1.0)
        if self.fine:
            self.gather_ise_fine()
        else:
            self.gather_ise()
        self.avepdfs = 0.0  # dummy
        self.fracise = 0.0  # dummy

    def printall(self, fine=False, lv=False, cond=False):
        self.fine = fine
        self.lv = lv
        self.cond = cond
        timelist = [str(t) + "m" for t in self.times]
        for tidx in range(1, 21):
            print(tidx)
            trydir = self.head + str(tidx) + "/"
            if self.fine and not self.lv:
                outfile = trydir + "ise_searched_rev"
                rfile = trydir + "result.txt_ise_rev"
            elif self.lv and not self.fine:
                outfile = trydir + "ise_searched_lv"
                rfile = trydir + "result.txt_ise_lv"
            elif self.cond and not self.lv and not self.fine:
                outfile = trydir + "ise_searched_withcond"
                rfile = trydir + "result.txt_ise_withcond"
            elif self.lv and self.fine:
                outfile = trydir + "ise_searched_rev_lv"
                rfile = trydir + "result.txt_ise_rev_lv"
            else:
                outfile = trydir + "ise_searched"
                rfile = trydir + "result.txt_ise"
            with open(outfile, 'w') as f, open(rfile, 'w') as r:
                for idx, tname in enumerate(timelist):
                    print(tname)
                    self.workdir = trydir + tname + "/"
                    if self.lv:
                        self.log = self.workdir + "std-ise-lv-" + str(tidx) + "_"\
                            + tname + ".log"
                    elif self.cond:
                        self.log = self.workdir + "std-ise-withcond-" + str(tidx) + "_"\
                            + tname + ".log"
                    else:
                        self.log = self.workdir + "std-ise-" + str(tidx) + "_"\
                            + tname + ".log"
                    self.reflog = self.workdir + "std-" + tname + ".log"
                    if self.shift is not None:
                        self._shift = self.shift[idx, :]
                    self.get_ise_from_one_dir()
                    print('{0}, {1}, {2}, {3}'.format(self.n, self.minise,
                                                      self.avepdfs,
                                                      self.fracise), file=f)
                    if idx == 0:
                        r.write("optimization results: n, optid_x, optid_y, " +
                                "optid_z, optid_w, 1/n, 1/(optidnx_x*optid_y" +
                                "*optid_z*optid_w)\n")
                    r.write("%e %d %d %d %d %e %e\n" %
                            (
                             self.n, self.minbw[0], self.minbw[1],
                             self.minbw[2], self.minbw[3],
                             1/(self.n*1.0), 1/(np.prod(self.minbw)*1.0)
                            )
                            )
                    if idx == len(timelist) - 1:
                        r.write("The end of data")


#def samplerun():
#    times = [2, 4, 6, 8, 10, 12, 14, 16, 24, 28, 32, 34, 67, 112, 168, 180, 284]
#    head = "/home/kazu/desktop/200701/orthotope_again_ddscs/try"
#    project = gather_ise_search(times, head)
#    project.printall()


#samplerun()
