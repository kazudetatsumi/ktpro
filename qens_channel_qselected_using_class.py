#!/usr/bin/env python
# example command: qens_qselected_using_class.py  6204 000025io 80 1 Boxcar
import sys
sys.path.append("/home/kazu/ktpro")
import qens_channel_class as qc


def run():
    print(sys.argv)
    if len(sys.argv) >= 2:
        runno = sys.argv[1]
        print("runno =", runno)
    else:
        print("using default runno 6204")
        runno = "6204"
    if len(sys.argv) >= 3:
        dirname = sys.argv[2]
        print("dirname =",dirname)
    else:
        print("using default dirname 000025io")
        dirname = "000025io"
    if len(sys.argv) >= 4:
        M = int(sys.argv[3])
        print("M =", M)
    else:
        print("using default M 80")
        M = 80
    if len(sys.argv) >= 5:
        winparam = float(sys.argv[4])
        print("winparam =", winparam)
    else:
        print("using default winparam 1")
        winparam = float(1)
    if len(sys.argv) >= 6:
        WinFunc = sys.argv[5]
        print("WinFunc =", WinFunc)
    else:
        print("using default WinFunc Boxcar")
        WinFunc = "Boxcar"

    datadir = "/home/kazu/desktop/210108/Tatsumi/srlz/"+dirname+"/"
    outdir = "/home/kazu/desktop/210108/Tatsumi/pickles/"+dirname+"/"
    save_file = datadir + "/run"+runno+"united_spectra.pkl"
    #output_file = outdir + "qens_run"+runno+"united_kde_results_on_data_qsel.pkl"
    output_file = "./qens_run"+runno+"united_kde_results_on_data_qsel.pkl"
    figname = "qens_out_"+runno+"_"+dirname+"_"+str(M)+"_"+str(winparam)+"_"+WinFunc+".png"
    proj = qc.qens(datadir, save_file, qsel=True, winparam=winparam, M=M,
                   WinFunc=WinFunc, figname=figname, showplot=False)
    proj.select_spectra()
    proj.add_shift()
    proj.run_ssvkernel()
    proj.save_output(output_file)


run()
