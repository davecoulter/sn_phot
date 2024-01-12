from matplotlib import pyplot as plt
import numpy as np
import sncosmo
from astropy.table import Table
import os,shutil,glob,sys
from models import BAYESNSource
from classify import *
import corner
from astropy.cosmology import LambdaCDM
import pdb
from stardust import classify
import pickle

cosmo = LambdaCDM(Om0=0.3, Ode0=0.7, H0=70.0)

do_classify = True

input_path = "./jades_phot"
output_path = "./fit_output"

# input_files = glob.glob("./jades_phot/*_sncosmo_phot.txt")
input_files = glob.glob("{input_path}/tr10_sncosmo_phot.txt".format(input_path=input_path))
z_file = Table.read("{input_path}/JADES_Sources1_Redshifts.txt".format(input_path=input_path), format='ascii')

for f in input_files:

    # book keeping
    file_name = os.path.basename(f)
    pickle_file = f.replace("txt", "pkl")
    pickle_path = "{output_path}/{pickle_file}".format(output_path=output_path, pickle_file=pickle_file)

    lc_plot_file = f.replace("txt", "png")
    lc_plot_path = "{output_path}/{lc_plot_file}".format(output_path=output_path, lc_plot_file=lc_plot_file)

    corner_plot_file = f.replace("txt", "corner.png")
    corner_plot_path = "{output_path}/{corner_plot_file}".format(output_path=output_path, corner_plot_file=corner_plot_file)


    lc = Table.read(f,format='ascii')
    lc = lc[~np.isnan(lc['flux'])]

    t0 = (np.min(lc['mjd'])-100, np.max(lc['mjd'])+100)

    # def classify(sn, zhost=1.491, zhosterr=0.003, t0_range=None,
    #              zminmax=[1.488,1.493], npoints=100, maxiter=10000,
    #              templateset='SNANA', excludetemplates=[],
    #              nsteps_pdf=101, priors={'Ia':0.24, 'II':0.57, 'Ibc':0.19},
    #              inflate_uncertainties=False,
    #              verbose=True):

    output = None 
    if do_classify:
        output = classify.classify(lc, zhost=3.6, zhosterr=0.001, t0_range=t0, zminmax=[3.599,3.601])        
        pickle.dump(output, open(pickle_path, "wb"))
    else:
        output = pickle.load(open(pickle_path, "rb"))

    model_name = output["bestmodel"]
    mod = output[model_name]["fit"]
    peak_mag = mod.source_peakabsmag("bessellb", "ab")
    print("\nBest model: %s; Peak mag: %0.3f [mag]" % (model_name, peak_mag))

    lc_fit = output[model_name]["sn"]
    sncosmo.plot_lc(lc_fit, mod)

    plt.savefig(lc_plot_path)
    plt.close()

    res = output[model_name]["res"]
    corner.corner(res.samples, weights=res.weights)
    plt.savefig(corner_plot_path)


