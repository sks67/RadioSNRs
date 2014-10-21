#-----------------------------------------------------------#
# This will be the most optimized, robust, efficient version
# of the driver program.
#-----------------------------------------------------------#

import numpy as np
import matplotlib.pyplot as py
import scipy.interpolate as interp
import scipy.integrate as integ
import trumck_diam
import ismdensgen as dgen
import time as tm
import datetime as dt
import os


a = np.array([0.5,0.75,1.0,1.25,1.5])
thick_array = np.array([100,150,200,250,300])*1.0
frac_array = np.array([0.1,0.3,0.5,0.7,0.9])
ratelmc = np.array([0.8,0.9,1.0,2.0,3.0])*1.0e-3*15.0
#a = np.array([1.4])
#thick_array = np.array([100])*1.0
#frac_array = np.array([0.3])
#ratelmc = np.array([3.0])*1.0e-3

likelihood = np.zeros(thick_array.size*frac_array.size*ratelmc.size*a.size)
delta_rate = 0.25
delta_a = 0.25
xx=0
start = tm.time()
for a_ind,a_elem in enumerate(a):
    print a_elem
                
    h1_1a = dgen.type1a_densities()
    h1_cc = dgen.ccsn_densities(a_elem)
    mass_cc,mass_1a = dgen.superMasses(h1_1a,h1_cc)
    ek = dgen.superEnergies(h1_1a,h1_cc)
    for j in range(thick_array.size):
        for i in range(frac_array.size):

        #start = tm.time()
            frac = frac_array[i]
            for rate_ind,rate_elem in enumerate(ratelmc):

                rate_1a_lmc = (frac/(1+frac))*rate_elem
                rate_cc_lmc = (1/(1+frac))*rate_elem
                time_1a = dgen.timegen(rate_1a_lmc)
                time_cc = dgen.timegen(rate_cc_lmc)

                #print 'h1_1a.size = ',h1_1a.size
                #print 'h1_cc.size = ',h1_cc.size
                #print 'time_1a.size = ',time_1a.size
                #print 'time_cc.size = ',time_cc.size
                time,h1,ejmas,energ,nprof = dgen.create_snarray(time_1a,time_cc,mass_1a,mass_cc,h1_1a,h1_cc,ek)
                likelihood[xx] = trumck_diam.radiolightcurve(thick_array[j],h1,time,ejmas,energ,nprof)
               # print a_elem,'\t',thick_array[j],'\t',frac,'\t',rate_elem,'\t',likelihood[xx]
               # print likelihood[xx]
                xx=xx+1
                
print 'RUNTIME = ',str(dt.timedelta(seconds=tm.time()-start))


userdoc = os.path.join(os.getcwd(),'DataAnalysis')                                                                                                      
np.savetxt(os.path.join(userdoc,'checksnaps_easier_likelihood.txt'),likelihood)
np.savetxt(os.path.join(userdoc,'checksnaps_easier_thickness.txt'),thick_array)
np.savetxt(os.path.join(userdoc,'checksnaps_easier_fraction.txt'),frac_array)
np.savetxt(os.path.join(userdoc,'checksnaps_easier_snrate.txt'),ratelmc)
