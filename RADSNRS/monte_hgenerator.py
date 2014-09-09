################################################################
#                                                              #
# PURPOSE: Generate neutral hydrogen densities where SNe expl- #
#          ode based on the LMC H1 distribution.               #
# AUTHOR: Sumit (duh!)
#                                                              #
################################################################

import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as fits
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
import scipy.interpolate as interp

a = np.genfromtxt('trumck_inputs.txt',dtype='str',comments='#',skip_header=4)
path_file = a[0]
temp_data = fits.getdata(path_file)

#Cleaning the FITs file                                                                                                                                                            
temp_data2 = temp_data[~np.isnan(temp_data)]
r_harray = temp_data2[temp_data2>0.]

#Creating the histogram
n,bins,patches = plt.hist(r_harray,6000,normed=True,histtype='step')
bin_vals = bins[1:]-bins[:-1] #Way to get the bin-sizes and check the normalization. woohoo!
print np.sum(bin_vals*n)
custm = stats.rv_discrete(name='custm',values = (bins[1:],n*bin_vals))
#Check if the cdf is working
plt.plot(bins,custm.cdf(bins))
plt.xscale('log')
plt.xlim(1.0e8,1.0e23)
plt.ylim(-0.5,1.5)
plt.show()

#Random number generation
num_rvs = 1000000
cdf_table = custm.cdf(bins)
prob_rand = np.random.random(num_rvs)
f = interp.interp1d(cdf_table,bins)
randh_new = f(prob_rand)

np.savetxt('trumck_ia_dens.txt',randh_new)
