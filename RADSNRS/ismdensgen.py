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
import scipy.integrate as integ
import ConfigParser
import random as rn


Config = ConfigParser.ConfigParser()
Config.read('config.ini')
path_to_file = Config.get('InputFiles','path')
hfile = Config.get('InputFiles','hydrogenmap')
temp_data = fits.getdata(path_to_file+hfile)
pc = 3.086e18 # in cm
lmc_thick = 100*pc #according to http://adsabs.harvard.edu/cgi-bin/bib_query?1973MNRAS.163..163W


#Cleans the nans and other contaminants in the data. PUT THIS IN DOCSTRINGS SOON!
def cleandata(temp_data):
    temp_data = temp_data[0] if temp_data.ndim > 2.0 else temp_data
    temp_data2 = temp_data[~np.isnan(temp_data)]
    return temp_data2[temp_data2>0.]

r_harray = cleandata(temp_data)
n,bins = np.histogram(r_harray,bins=6000,normed=True)
#Getting the discrete pdf from combining Schmidt Kennicutt Law and LMC's HI distribution

def ccsn_densities(a):
    nh = r_harray/lmc_thick
    n_h_min = np.amin(nh)
    n_h_max = np.amax(nh)
    n_h = np.linspace(n_h_min,n_h_max,6000)
    norm = (a+1.0)/(n_h_max**(a+1.0) - n_h_min**(a+1.0))
    pdf_nh = norm*(n_h**a)*n
    norm_ccsn = 1.0/np.sum(pdf_nh*n_h)
    pdf_nh_act = norm_ccsn*pdf_nh
#Normalizing the discrete pdf and generating random densities
    custm = stats.rv_discrete(name='custm',values=(n_h,pdf_nh_act*n_h))
    num_rvs = 5000000
    cdf_table = custm.cdf(n_h)
    prob_rand = np.random.random(num_rvs)
    f = interp.interp1d(cdf_table,n_h)
    nh_new = f(prob_rand)
    return nh_new*lmc_thick

def type1a_densities():
    bin_vals = bins[1:] - bins[:-1]
    custm = stats.rv_discrete(name = 'custm', values = (bins[1:], n*bin_vals))
    num_rvs = 5000000
    cdf_table = custm.cdf(bins)
    prob_rand = np.random.random(num_rvs)
    f = interp.interp1d(cdf_table, bins)
    randh_new = f(prob_rand)
    return randh_new

def densityplots(snia,sncc):
    binarray = np.logspace(18,23,2000)
    plt.hist(snia,bins=binarray,histtype='step')
    plt.hist(sncc,bins=binarray,color='r',histtype='step')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(1,1.0e5)
    plt.xlim(1.0e18,1.0e23)
    plt.show()

#-----------------------------------------------------------#                                                                                                                    
#     SUPERNOVA EJECTA MASSES AND KINETIC ENERGIES          #                                                                                                                    
#-----------------------------------------------------------#                                                                                                                   

    #CCSN Masses picked based on high-end IMF.                                                                                                          
def superMasses(h1_1a,h1_cc):
    mass_cc = np.ones(h1_cc.size)*12.0
    mass_1a = np.ones(h1_1a.size)*1.4
    return (mass_cc,mass_1a)

def superEnergies(h1_1a,h1_cc):
    rv = np.random.normal(loc = 51.0, scale = 0.3, size = (h1_1a.size+h1_cc.size)) #0.4 variance ensures Hypernova, PISNs etc occur 1/1000 times "normal" CCSN (Janka 2012)     
    return 10**(rv-51.0)

#-----------------------------------------------------------#                                                                                         
#   MERGING SUBROUTINE FOR BIRTH TIMES AND H-DENSITY        #                                                                                          
#-----------------------------------------------------------#                                                                                         

def create_snarray(t_1a, t_cc,mass_1a,mass_cc,h1_1a,h1_cc,ek):
    t = np.concatenate([t_1a, t_cc])
    m = np.concatenate([mass_1a[0:t_1a.size],mass_cc[0:t_cc.size]])
    h = np.concatenate([h1_1a[0:t_1a.size],h1_cc[0:t_cc.size]])
    n_1a = np.ones(t_1a.size)*7.0
    n_cc = np.ones(t_cc.size)*12.0
    n = np.concatenate([n_1a, n_cc])

    args = np.argsort(t)
    t = t[args]
    h = h[args]
    m = m[args]
    sn_ek = ek[args]
    n = n[args]
    return (t, h, m, ek, n)

#-----------------------------------------------------------#                                                                                 
#   POISSON PROCESS BIRTH TIME GENERATOR SUBROUTINE         #                                                                                            
#-----------------------------------------------------------#                                                                                           

def timegen(snrate):
   t = np.arange(1.0e7)
   prob = np.random.rand(t.size)
   return t[np.where(prob<=snrate)]
   


#np.savetxt('trumck_cc_dens.txt',ccsn_densities())
#np.savetxt('trumck_ia_dens.txt',type1a_densities())
