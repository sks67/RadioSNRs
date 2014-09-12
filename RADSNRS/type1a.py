#-----------------------------------------------------------#
# This will be the most optimized, robust, efficient version
# of the driver program.
#-----------------------------------------------------------#

import numpy as np
import matplotlib.pyplot as py
import scipy.interpolate as interp
import scipy.integrate as integ
import trumck
import ismdensgen as dgen
import time as tm
import datetime as dt
import os


a = np.array([1.25,1.5])
thick_array = np.array([50])
frac_array = np.array([0.3])
ratelmc = np.array([3.0e-2,4.0e-2])
delta_rate = 0.25
#delta_a = a[1]-a[0]
for j in range(thick_array.size):
    for i in range(frac_array.size):

        xx=0
        likelihood = np.zeros(ratelmc.size*a.size)
        #start = tm.time()
        frac = frac_array[i]
        for rate_ind,rate_elem in enumerate(ratelmc):

            rate_1a_lmc = (frac/(1+frac))*rate_elem
            rate_cc_lmc = (1/(1+frac))*rate_elem
            time_1a = dgen.timegen(rate_1a_lmc)
            time_cc = dgen.timegen(rate_cc_lmc)
            for a_ind,a_elem in enumerate(a):
                
                h1_1a = dgen.type1a_densities()
                h1_cc = dgen.ccsn_densities(a_elem)
                mass_cc,mass_1a = dgen.superMasses(h1_1a,h1_cc)
                ek = dgen.superEnergies(h1_1a,h1_cc)
                time,h1,ejmas,energ,nprof = dgen.create_snarray(time_1a,time_cc,mass_1a,mass_cc,h1_1a,h1_cc,ek)
                likelihood[xx] = trumck.radiolightcurve(thick_array[j],h1,time,ejmas,energ,nprof)
                print a_elem,'\t',rate_elem,'\t',likelihood[xx]
                xx=xx+1

        print np.sum(likelihood)
        

#userdoc = os.path.join(os.getcwd(),'DataAnalysis')                                                                                                      
#np.savetxt(os.path.join(userdoc,'checksnaps_easier_likelihood.txt'),likelihood)
#np.savetxt(os.path.join(userdoc,'checksnaps_easier_thickness.txt'),thick_array)
#np.savetxt(os.path.join(userdoc,'checksnaps_easier_fraction.txt'),frac_array)
#np.savetxt(os.path.join(userdoc,'checksnaps_easier_snrate.txt'),ratelmc)
