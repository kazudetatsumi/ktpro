#!/usr/bin/env python
# gather lowest ise in try directries.
# Kazuyoshi TATSUMI 2020.
import numpy as np
import h5py
import subprocess


class gather_ise_search:
    def __init__(self, times, head):
        self.times = times
        self.head = head

    def gather_ise(self):
        line = subprocess.getoutput("tail --line=1 " + self.logfile)
        print(line, self.logfile)
        values = line.split("(")
        values2 = values[1].split(")")
        self.binwidths = np.array(values2[0].split(","), dtype='int32')
        self.ise = float(values2[1].split(" ")[5])

    def read_h5py(self, outfile, term):
        f = h5py.File(outfile, 'r')
        return np.array(f[term])

    def get_ise_from_one_dir(self):
        self.n = np.sum(self.read_h5py(self.workdir+"/eliminated_data.hdf5",
                                       "data4")*1.0)
        self.gather_ise()
        self.avepdfs = 0.0  # dummy
        self.fracise = 0.0  # dummy

    def printall(self):
        timelist = [str(t) + "m" for t in self.times]
        for tidx in range(1, 21):
            trydir = self.head + str(tidx) + "/"
            outfile = trydir + "ise_searched"
            rfile = trydir + "result.txt_ise"
            with open(outfile, 'w') as f, open(rfile, 'w') as r:
                for idx, tname in enumerate(timelist):
                    self.workdir = trydir + tname + "/"
                    self.logfile = self.workdir + "std-ise-" + str(tidx) + "_"\
                        + tname + ".log"
                    self.get_ise_from_one_dir()
                    print('{0}, {1}, {2}, {3}'.format(self.n, self.ise,
                                                      self.avepdfs,
                                                      self.fracise), file=f)
                    if idx == 0:
                        r.write("optimization results: n, optid_x, optid_y, " +
                                "optid_z, optid_w, 1/n, 1/(optidnx_x*optid_y" +
                                "*optid_z*optid_w)\n")
                    r.write("%e %d %d %d %d %e %e\n" %
                            (
                             self.n, self.binwidths[0], self.binwidths[1],
                             self.binwidths[2], self.binwidths[3],
                             1/(self.n*1.0), 1/(np.prod(self.binwidths)*1.0)
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
