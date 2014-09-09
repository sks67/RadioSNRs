#-----------------------------------------------------------#
# This will be the most optimized, robust, efficient version
# of the driver program.
#-----------------------------------------------------------#

import numpy as np
import matplotlib.pyplot as py
import scipy.interpolate as interp
import scipy.integrate as integ
import trumck
import time as tm
import datetime as dt
import random as rn
import os
import ConfigParser


Config = ConfigParser.ConfigParser()
Config.read('config.ini')
nlmc = Config.getint('ObsParameters','numberofsnrs')


#LATER....create h1_{} arrays of the size of the largest time array i.e. for rate[-1]
#and then use subsets for smaller snrates

h1_1a = np.loadtxt('trumck_ia_dens.txt')
h1_cc = np.loadtxt('trumck_cc_dens.txt')

#-----------------------------------------------------------#
#     SUPERNOVA EJECTA MASSES AND KINETIC ENERGIES          #
#-----------------------------------------------------------#

#CCSN Masses picked based on high-end IMF.
xmin = 2 #solar masses
xmax = 30 #solar masses
n=2.35
norm = (1.0-n)/(xmax**(1.0-n) - xmin**(1.0-n))
pdf_x = lambda x: norm*x**(-n)

#Check Normalization
y,err = integ.quad(pdf_x,xmin,xmax)
cdf_x = lambda x: integ.quad(pdf_x,xmin,x)

#Random Variables by interpolation
x = np.linspace(xmin,xmax,5000)
cdf_table = np.zeros_like(x)
for i in range(x.size):
    y,err = cdf_x(x[i])
    cdf_table[i]=y

f = interp.interp1d(cdf_table,x)
cdfnew = np.random.random(size=h1_cc.size)
mass_cc = f(cdfnew)

#For Type Ia, Scalzo et al 2014 (Nearby SN Factory)
mass_1a = np.random.uniform(0.9,1.4,size=h1_1a.size)

rv = np.random.normal(loc=51.0,scale=0.28,size=(h1_1a.size+h1_cc.size)) 
#0.4 variance ensures Hypernova, PISNs etc occur 1/1000 times "normal" CCSN (Janka 2012)

ek = 10**(rv-51.0)

#-----------------------------------------------------------#
#   MERGING SUBROUTINE FOR BIRTH TIMES AND H-DENSITY        #
#-----------------------------------------------------------#

def create_snarray(t_1a,t_cc):
    t = np.concatenate([t_1a,t_cc])
    m = np.concatenate([mass_1a[0:t_1a.size],mass_cc[0:t_cc.size]])
    h = np.concatenate([h1_1a[0:t_1a.size],h1_cc[0:t_cc.size]])
    n_1a = np.ones(t_1a.size)*7.0
    n_cc = np.ones(t_cc.size)*12.0
    n = np.concatenate([n_1a,n_cc])

    args = np.argsort(t)
    t = t[args]
    h = h[args]
    m = m[args]
    sn_ek = ek[args]
    n = n[args]
    return (t,h,m,ek,n)

#-----------------------------------------------------------#
#   POISSON PROCESS BIRTH TIME GENERATOR SUBROUTINE         #
#-----------------------------------------------------------#

tmax = 4.0e7

def timegen(snrate):
   time = []
   j = 0
   count = 0
   for i in xrange(int(tmax)):
       prob = rn.random()
       if prob > snrate:
           continue
       else:
           count=count+1
           time.append(i)
           j=j+1

   return np.array(time)

#-----------------------------------------------------------#
#   EXPLORING SCALE HEIGHT-FRACTION PARAMETER SPACE         #
#-----------------------------------------------------------#


thick_array = np.array([100.0])
frac_array = np.array([0.8])
ratelmc = np.array([4.0e-3])
xx=0
likelihood = np.zeros(thick_array.size*frac_array.size)
for j in range(thick_array.size):
    for i in range(frac_array.size):
        start = tm.time()
        frac = frac_array[i]

        rate_1a_lmc = (frac/(1+frac))*ratelmc[0]
        rate_cc_lmc = (1/(1+frac))*ratelmc[0]
        
        time_1a = timegen(rate_1a_lmc)
        time_cc = timegen(rate_cc_lmc)
   
        time,h1,ejmas,energ,nprof = create_snarray(time_1a,time_cc)
        
        likelihood[xx] = trumck.radiolightcurve(thick_array[j],h1,time,ejmas,energ,nprof)
        
        print likelihood[xx],'\t','RUNTIME: ',dt.timedelta(seconds=tm.time()-start)
        xx=xx+1

#userdoc = os.path.join(os.getcwd(),'DataAnalysis')                                                                                                      
#np.savetxt(os.path.join(userdoc,'checksnaps_easier_likelihood.txt'),likelihood)
#np.savetxt(os.path.join(userdoc,'checksnaps_easier_thickness.txt'),thick_array)
#np.savetxt(os.path.join(userdoc,'checksnaps_easier_fraction.txt'),frac_array)
#np.savetxt(os.path.join(userdoc,'checksnaps_easier_snrate.txt'),ratelmc)
