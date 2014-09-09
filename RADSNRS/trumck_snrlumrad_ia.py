#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Testing the light curve model for a single supernova remnant in a given#
# ISM density.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

import datetime as dt
import numpy as np
import matplotlib.pyplot as py
import math as mt
import scipy.interpolate as interp
import subprocess
import sys

#random locations of densities and times for snrate = 2.0e3
lmcthick = 3.086e20 #gas thickness according to MB10.

#mass = np.array([1.4,5.0,10.0,15.0])
lum = np.zeros(800) #800 is the size of lluni2 for the sedov taylor phases
tim = np.zeros(800) #800 is the size of lluni2 for the sedov taylor phases

rad = np.zeros(800) #I added this to check radius vs time for a single SN. 


nt2 = 800
lluni2 = np.zeros(nt2)
mu=1.4
n0 = 1.0 #dimensionless, multiplied by 1cm^3
mp = 1.67e-24
ismrho = n0*mu*mp #in units of cm^-3
e51 = 2.0 #in units of 10^51 ergs of energy released per sne
mej = 1.4 #in units of solar masses

#Characteristic scales
tch = 423*(e51**(-0.5))*(mej**(5.0/6.0))*(n0**(-1.0/3.0)) #years
rch = 3.07*(mej**(1.0/3.0))*(n0**(-1.0/3.0)) #pcs
vch = 7090*(e51**0.5)*(mej**(-0.5)) #km/s

#variables 
tstar_st = 0.732 #beggining of dimensionless sedov taylor phase
t_st0 = tstar_st*tch


#________________________[ n=7 ]______________________________# 

#--------------------EJECTA PHASE-------------------------#
t_ed = np.logspace(0.1,np.log10(t_st0),1000)
tstar_ed = t_ed/tch
vstar_ed = 0.606*tstar_ed**(-(3.0/7.0))
rstar_ed = 1.06*tstar_ed**(4.0/7.0)
v_ed = vstar_ed*vch  #in km/s
r_ed = rstar_ed*rch  #in parsecs
#--------------------------------------------------------#

#--------------------SEDOV PHASE-------------------------#
t_st = np.arange(nt2)/200.0 + np.log10(t_st0)
tstar = (10**t_st)/tch
vstar = 0.569*((1.42*tstar - 0.312)**(-3.0/5.0))
rstar = (1.42*tstar - 0.312)**(2.0/5.0)
v_st = vstar*vch  #in km/s
r_st = rstar*rch  #in parsecs
#--------------------------------------------------------#


epsb = 0.01
alpha = 1.0 
epse = epsb*alpha
pp = 3.0
compf = 4.0
ff = 0.38
rho0 = ismrho/((1.0e-24)*mu*1.67)
nei = 1.14
cl = 6.27e18
c5 = 7.52e-24
c6 = 7.97e-41
dist = 50*1000*3.086e18
nu = 1.4e9

#-------------------------------ED STAGE LUMINOSITY AND RADIUS---------------------------------#
bism = 9e-6*rho0**0.47    #magnetic field for galaxies scales as 0.47 (Krutcher et al. 1999)
rr1 = r_ed*(3.086e18) #LAURA - cm
vv1 = v_ed*1.0e5 #LAURA - cm/s
bb0 = np.sqrt(8.0*mt.pi*epsb*ismrho*vv1*vv1)
bb1 = np.sqrt(8.0*mt.pi*epsb*ismrho*vv1*vv1)
goat1 = np.where(bb1<(4.0*bism))
ngoat1 = goat1[0].size
if ngoat1>0.0:
    bb1[goat1]=4.0*bism
    
me = 9.11e-28  #in grams
c = 3.0e10
gammam = (mu*epse*(pp-2)*((vv1/c)**2)*mp)/((pp-1)*compf*nei*me)
goat = np.where(gammam<1.0)
ngoat = goat[0].size

if ngoat>0.0:
    gammam[goat]=1.0
    
elow = gammam*me*c**2  #Eq. 10 Chevalier, 98
n_0 = (alpha*(bb0**2)*(pp-2)*(elow**(pp-2.0)))/(8.0*mt.pi)
ss11 = (4.0/3.0)*ff*rr1
nu11 = 2.0*cl*((ss11*c6*n_0)**(2.0/(pp+4.0)))*(bb1**((pp+2.0)/(pp+4.0)))
ss1 = (c5/c6)*(bb1**(-0.5))*((nu11/(2.0*cl))**2.5)
jj1 = ((nu/nu11)**2.5)*(1.0-np.exp(-1.0*(nu/nu11)**(-1.0*(pp+4.0)/2.0)))
f_nu1 = (ss1*jj1*mt.pi*rr1**2)/(1.0e-26*dist**2)
lluni = f_nu1*1.0e-26*4.0*mt.pi*dist**2
#----------------------------------------------------------------------------------------------#


#-------------------------------ST STAGE LUMINOSITY AND RADIUS---------------------------------#

bism = 9e-6*rho0**0.47    #magnetic field for galaxies scales as 0.47 (Krutcher et al. 1999)
rr = r_st*(3.086e18) #LAURA - cm
vv = v_st*1.0e5 #LAURA - cm/s
bb0 = np.sqrt(8.0*mt.pi*epsb*ismrho*vv*vv)
bb = np.sqrt(8.0*mt.pi*epsb*ismrho*vv*vv)
goat = np.where(bb<(4.0*bism))
ngoat = goat[0].size
if ngoat>0.0:
    bb[goat]=4.0*bism
    
me = 9.11e-28  #in grams
c = 3.0e10
gammam = (mu*epse*(pp-2)*((vv/c)**2)*mp)/((pp-1)*compf*nei*me)
goat = np.where(gammam<1.0)
ngoat = goat[0].size


if ngoat>0.0:
    gammam[goat]=1.0
    
elow = gammam*me*c**2  #Eq. 10 Chevalier, 98
n_0 = (alpha*(bb0**2)*(pp-2)*(elow**(pp-2.0)))/(8.0*mt.pi)
ss1 = (4.0/3.0)*ff*rr
nu1 = 2.0*cl*((ss1*c6*n_0)**(2.0/(pp+4.0)))*(bb**((pp+2.0)/(pp+4.0)))
ss = (c5/c6)*(bb**(-0.5))*((nu1/(2.0*cl))**2.5)
jj = ((nu/nu1)**2.5)*(1.0-np.exp(-1.0*(nu/nu1)**(-1.0*(pp+4.0)/2.0)))
f_nu = (ss*jj*mt.pi*rr**2)/(1.0e-26*dist**2)
lluni2 = f_nu*1.0e-26*4.0*mt.pi*dist**2

#End of ST phase.        
vrad = 188*((e51/(n0*n0))**(0.07))
goat = np.where(v_st<=vrad)
lluni2[goat]=0.0
print 'v_rad = ',vrad,' km/s'
print 'numerical trad = ',(10**t_st)[goat[0]][0],' years'
print 'theoretical trad = ',(49300*(e51**0.22)*(n0**(-0.55))),' years'

#----------------------------------------------------------------------------------------------#

lum=np.concatenate([lluni,lluni2])
tim=np.concatenate([t_ed,10**t_st])
rad=np.concatenate([r_ed,r_st])


np.savetxt('trumck_snrlumrad_ia_tim.txt',tim)
np.savetxt('trumck_snrlumrad_ia_lum.txt',lum)
np.savetxt('trumck_snrlumrad_ia_rad.txt',rad)
np.savetxt('trumck_snrlumrad_ia_ted.txt',t_ed)
np.savetxt('trumck_snrlumrad_ia_ved.txt',v_ed)
np.savetxt('trumck_snrlumrad_ia_tst.txt',10**t_st)
np.savetxt('trumck_snrlumrad_ia_vst.txt',v_st)



fig = py.figure(1)
ax = fig.add_subplot(111)
ax.plot(tim,lum,'r-',lw=1.5)
ax.axvline(t_st0,linestyle='--',lw=1.5)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r't [years]',fontsize='18')
ax.set_ylabel(r'1.4 GHz Luminosity (ergs/s/Hz)',fontsize='18')

ax.set_xlim(100,1.0e6)
ax.set_ylim(1.0e19,1.0e26)
#ax.axis([10.0,1.0e5,1.0e18,1.0e24])
ax.legend()
np.savetxt('snrlum_e51_05.txt',lum)
np.savetxt('snrtim_e51_05.txt',tim)



fig2 = py.figure(2)
ax2 = fig2.add_subplot(111)
ax2.axvline(t_st0,linestyle='--',lw=1.5)
ax2.plot(tim,rad,'r-',lw=1.5)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel(r't [years]',fontsize='18')
ax2.set_ylabel(r'Radius [pc]',fontsize='18')
ax2.set_xlim(100,1.0e6)

ax2.legend(loc=4)

fig3 = py.figure(3)
ax3 = fig3.add_subplot(111)
ax3.axvline(t_st0,linestyle='--',lw=1.5)
ax3.plot(t_ed,v_ed,'r-',10**t_st,v_st,'r-',lw=1.5)
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xlabel(r't [years]',fontsize='18')
ax3.set_ylabel(r'Velocity [km/s]',fontsize='18')
ax3.set_xlim(100,1.0e6)

ax3.legend(loc=1)

py.tick_params(labelsize='16')
py.tight_layout()
py.show()
