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

input_files = glob.glob("./JADES_fits/*_sncosmo_phot.txt")
z_file = Table.read("./JADES_fits/JADES_Sources1_Redshifts.txt", format='ascii')

for f in input_files:
    lc = Table.read("./JADES_Sources1/tr10_sncosmo_phot.txt",format='ascii')
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
    	pickle.dump(output, open("./JADES_Sources1/tr10_sncosmo_phot.pkl", "wb"))
    else:
    	output = pickle.load(open("./JADES_Sources1/tr10_sncosmo_phot.pkl", "rb"))


model_name = output["bestmodel"]
mod = output[model_name]["fit"]
peak_mag = mod.source_peakabsmag("bessellb", "ab")
print("Best model: %s; Peak mag: %0.3f [mag]" % (model_name, peak_mag))

lc_fit = output[model_name]["sn"]
sncosmo.plot_lc(lc_fit, mod)

plt.savefig("./JADES_Sources1/tr10_snii_best_fit.png")
plt.close()

try:
	res = output[model_name]["res"]
	corner.corner(res.samples, weights=res.weights)
	plt.savefig("./JADES_Sources1/tr10_snii_best_fit_corner.png")
except:
	pdb.set_trace()
	test = 1

