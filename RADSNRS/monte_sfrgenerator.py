##################################################################
# 
#  PURPOSE: Generates HI based on Schmidt Kennicutt Law, but with
#          a slightly different technique, and probably more corr-
#          ect.
#
##################################################################

import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as fits
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
import scipy.interpolate as interp

import ConfigParser
Config = ConfigParser.ConfigParser()
Config.read('config.ini')
path_to_file = Config.get('InputFiles','path')
hfile = Config.get('InputFiles','hydrogenmap')
lmch1_file = Config.get('InputFiles','density')

temp_data = fits.getdata(path_to_file+hfile)
lmc_h1 = np.loadtxt(path_to_file+lmch1_file)

#Cleans the nans and other contaminants in the data. PUT THIS IN DOCSTRINGS SOON!
def cleandata(temp_data):
    temp_data = temp_data[0] if temp_data.ndim > 2.0 else temp_data
    temp_data2 = temp_data[~np.isnan(temp_data)]
    return temp_data2[temp_data2>0.]

a = 1.0
pc = 3.086e18 # in cm
lmc_thick = 100*pc #according to http://adsabs.harvard.edu/cgi-bin/bib_query?1973MNRAS.163..163W


r_harray = cleandata(temp_data)
n,bins = np.histogram(r_harray,bins=6000,normed=True)

#Getting the discrete pdf from combining Schmidt Kennicutt Law and LMC's HI distribution

def ccsn_densities():
    nh = r_harray/lmc_thick
    n_h_min = np.amin(nh)
    n_h_max = np.amax(nh)
    n_h = np.linspace(n_h_min,n_h_max,6000)
    norm = (a+1.0)/(n_h_max**(a+1.0) - n_h_min**(a+1.0))
    pdf_nh = norm*(n_h**a)*n
    norm_ccsn = 1.0/np.sum(pdf_nh*n_h)
    pdf_nh_act = norm_ccsn*pdf_nh
    print 'pdf sum = ', np.sum(pdf_nh_act*n_h)

#Normalizing the discrete pdf and generating random densities
    print np.sum(pdf_nh_act*n_h)
    custm = stats.rv_discrete(name='custm',values=(n_h,pdf_nh_act*n_h))
    num_rvs = 1000000
    cdf_table = custm.cdf(n_h)
    prob_rand = np.random.random(num_rvs)
    f = interp.interp1d(cdf_table,n_h)
    nh_new = f(prob_rand)
    return nh_new*lmc_thick

def type1a_densities():
    bin_vals = bins[1:] - bins[:-1]
    custm = stats.rv_discrete(name = 'custm', values = (bins[1:], n*bin_vals))
    num_rvs = 1000000
    cdf_table = custm.cdf(bins)
    prob_rand = np.random.random(num_rvs)
    f = interp.interp1d(cdf_table, bins)
    randh_new = f(prob_rand)
    return randh_new

np.savetxt('trumck_cc_dens.txt',ccsn_densities())
np.savetxt('trumck_ia_dens.txt',type1a_densities())
