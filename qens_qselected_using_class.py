#!/usr/bin/env python
import sys
sys.path.append("/home/kazu/ktpro")
import qens_class as qc


def run():
    datadir = "/home/kazu/desktop/210108/Tatsumi/srlz/000025io/"
    outdir = "/home/kazu/desktop/210108/Tatsumi/pickles/000025io/"
    save_file = datadir + "/run6204united_spectra.pkl"
    #sfile = datadir + "/run6204united_sspectra.pkl"
    output_file = outdir + "qens_run6204united_kde_results_on_data_qsel.pkl"
    #outputs_file = outdir + "qens_run6204united_kde_results_on_sdata_qsel.pkl"
    #
    #proj = qc.qens(datadir, save_file, sfile=sfile, qsel=True, optsm=True)
    proj = qc.qens(datadir, save_file, qsel=True)
    #proj.selected_spectra = proj.dataset['spectra']
    #proj.selected_energy = proj.dataset['energy']
    #print(proj.selected_energy[np.argmax(proj.selected_spectra)])
    #proj.de = proj.selected_energy[1] - proj.selected_energy[0]
    #proj.sde = proj.senergy[1] - proj.senergy[0]
    proj.select_spectra()
    #proj.get_xvec()
    proj.add_shift()
    proj.run_ssvkernel()
    proj.save_output(output_file)
    #proj.save_outputs(outputs_file)


run()
