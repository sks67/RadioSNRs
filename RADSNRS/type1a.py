#-----------------------------------------------------------#
# This will be the most optimized, robust, efficient version
# of the driver program.
#-----------------------------------------------------------#

import numpy as np
import matplotlib.pyplot as py
import scipy.interpolate as interp
import scipy.integrate as integ
import trumck2
import ismdensgen as dgen
import time as tm
import datetime as dt
import os



a = np.array([1.3])
#thick_array = np.array([100])*1.0
frac = 0.33
#ratelmc = np.array([3.0])*1.0e-3
#epsb = np.array([0.001])
thick_array = np.linspace(10,900,10)*1.0#np.array([50,250,450,650,850])*1.0
ratelmc = np.linspace(1.0,5.0,10)*1.0e-3*15.0#np.array([1.0,1.5,2.0,2.5,3.0])*1.0e-3*15.0
epsb = np.logspace(-3.0,-1.3,10)

likhood_lum = np.zeros(thick_array.size*epsb.size*ratelmc.size)
likhood_diam = np.zeros_like(likhood_lum)
likhood_dens = np.zeros_like(likhood_lum)
delta_rate = 0.25
delta_a = 0.25
xx=0
start = tm.time()

h1_1a = dgen.type1a_densities()
h1_cc = dgen.ccsn_densities(a)
mass_cc,mass_1a = dgen.superMasses(h1_1a,h1_cc)
ek = dgen.superEnergies(h1_1a,h1_cc)


for j in range(thick_array.size):
    print j
    for rate_ind,rate_elem in enumerate(ratelmc):
        rate_1a_lmc,rate_cc_lmc = (frac/(1+frac))*rate_elem,(1/(1+frac))*rate_elem
        time_1a,time_cc = dgen.timegen(rate_1a_lmc),dgen.timegen(rate_cc_lmc)
        for eb_i,eb in enumerate(epsb):
            time,h1,ejmas,energ,nprof = dgen.create_snarray(time_1a,time_cc,mass_1a,mass_cc,h1_1a,h1_cc,ek)
            likhood_lum[xx],likhood_diam[xx],likhood_dens[xx] = trumck2.radiolightcurve(thick_array[j],eb,h1,time,ejmas,energ,nprof)
            xx=xx+1

print 'RUNTIME = ',str(dt.timedelta(seconds=tm.time()-start))
        
        
userdoc = os.path.join(os.getcwd(),'DataAnalysis')                                                                                                      
np.savetxt(os.path.join(userdoc,'epsb_likhood_lum.txt'),likhood_lum)
np.savetxt(os.path.join(userdoc,'epsb_likhood_diam.txt'),likhood_diam)
np.savetxt(os.path.join(userdoc,'epsb_likhood_dens.txt'),likhood_dens)

#np.savetxt(os.path.join(userdoc,'checksnaps_easier_thickness.txt'),thick_array)
#np.savetxt(os.path.join(userdoc,'checksnaps_easier_fraction.txt'),frac_array)
#np.savetxt(os.path.join(userdoc,'checksnaps_easier_snrate.txt'),ratelmc)
