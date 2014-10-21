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
from mpl_toolkits.axes_grid1 import Grid

#random locations of densities and times for snrate = 2.0e3
lmcthick = 3.086e20 #gas thickness according to MB10.
n0_array = np.array([0.1,1.0,10.0])
e51_array = np.array([0.5,1.0,2.0])
mej_array = np.array([2.0,10.0,20.0])
initmatrix = np.vstack([n0_array,e51_array,mej_array])
colors = np.array(['r','k','b'])
#Create gridplot framework
fig = py.figure(1,figsize=(10,10))
grid = Grid(fig,rect=111,nrows_ncols=(3,3),axes_pad=0.0,label_mode='L')

def textstring(param,value):
    if param==0:
        return ('$n_0$','$cm^{-3}$')
    elif param==1:
        return ('$E_{51}$',' ')
    else:
        return ('$M_{ej}$','$M_{\odot}$')
        
for x,arr in enumerate(initmatrix):
    for ind,elem in enumerate(arr):
       
#mass = np.array([1.4,5.0,10.0,15.0])
            lum = np.zeros(800) #800 is the size of lluni2 for the sedov taylor phases
            tim = np.zeros(800) #800 is the size of lluni2 for the sedov taylor phases
            
            rad = np.zeros(800) #I added this to check radius vs time for a single SN. 
    
    
            nt2 = 800
            lluni2 = np.zeros(nt2)
            mu=1.4
            nh = 1.0
    #n0 = 1.0 #dimensionless, multiplied by 1cm^3
            mp = 1.67e-24
            n0 = elem if x==0 else initmatrix[0,1]
            ismrho = n0*mu*mp #in units of cm^-3
            e51= elem if x==1 else initmatrix[1,1]
            mej= elem if x==2 else initmatrix[2,1]
            

#Characteristic scales
            tch = 423*(e51**(-0.5))*(mej**(5.0/6.0))*(n0**(-1.0/3.0)) #years
            rch = 3.07*(mej**(1.0/3.0))*(n0**(-1.0/3.0)) #pcs
            vch = 7090*(e51**0.5)*(mej**(-0.5)) #km/s


#--------------------EJECTA PHASE-------------------------#                                                                                                                     
            tstar_st = 0.424
            t_st0 = tstar_st*tch
            t_ed = np.logspace(0.1,np.log10(t_st0),1000)
            tstar_ed = t_ed/tch
            vstar_ed = 0.593*tstar_ed**(-(3.0/7.0))
            rstar_ed = 1.038*tstar_ed**(4.0/7.0)
            v_ed = vstar_ed*vch  #in km/s                                                                                                                                             
            r_ed = rstar_ed*rch  #in par                                                                                                                                            
#--------------------SEDOV PHASE-------------------------#                                                                                                                          
            t_st = np.arange(nt2)/200.0 + np.log10(t_st0)
            tstar = (10**t_st)/tch
            vstar = 0.434*((1.42*tstar - 0.28)**(-3.0/5.0))
            rstar = (1.42*tstar - 0.28)**(2.0/5.0)
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
            print gammam
            goat = np.where(gammam<1.0)
            ngoat = goat[0].size
    
            if ngoat>0.0:
                gammam[goat]=1.0
        
            elow = gammam*me*c**2  #Eq. 10 Chevalier, 98
            n0 = (alpha*(bb0**2)*(pp-2)*(elow**(pp-2.0)))/(8.0*mt.pi)
            ss11 = (4.0/3.0)*ff*rr1
            nu11 = 2.0*cl*((ss11*c6*n0)**(2.0/(pp+4.0)))*(bb1**((pp+2.0)/(pp+4.0)))
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
            print gammam
            goat = np.where(gammam<1.0)
            ngoat = goat[0].size
    

            if ngoat>0.0:
                gammam[goat]=1.0
        
            elow = gammam*me*c**2  #Eq. 10 Chevalier, 98
            n0 = (alpha*(bb0**2)*(pp-2)*(elow**(pp-2.0)))/(8.0*mt.pi)
            ss1 = (4.0/3.0)*ff*rr
            nu1 = 2.0*cl*((ss1*c6*n0)**(2.0/(pp+4.0)))*(bb**((pp+2.0)/(pp+4.0)))
            ss = (c5/c6)*(bb**(-0.5))*((nu1/(2.0*cl))**2.5)
            jj = ((nu/nu1)**2.5)*(1.0-np.exp(-1.0*(nu/nu1)**(-1.0*(pp+4.0)/2.0)))
            f_nu = (ss*jj*mt.pi*rr**2)/(1.0e-26*dist**2)
            lluni2 = f_nu*1.0e-26*4.0*mt.pi*dist**2
    
#End of ST phase.        
            vrad = 200.
            goat = np.where(v_st<=vrad)
            lluni2[goat]=0.0
#----------------------------------------------------------------------------------------------#

          #  print lluni[-1],lluni2[0]
            lum=np.concatenate([lluni,lluni2])
          #  print t_ed[-1],10**t_st[0]
            tim=np.concatenate([t_ed,10**t_st])
          #  print r_ed[-1],r_st[0]
            rad=np.concatenate([r_ed,r_st])

  #  np.savetxt('trumck_snrlumrad_cc_tim.txt',tim)
  #  np.savetxt('trumck_snrlumrad_cc_lum.txt',lum)
  #  np.savetxt('trumck_snrlumrad_cc_rad.txt',rad)
   # np.savetxt('trumck_snrlumrad_cc_ted.txt',t_ed)
   # np.savetxt('trumck_snrlumrad_cc_ved.txt',v_ed)
    #np.savetxt('trumck_snrlumrad_cc_tst.txt',10**t_st)
    #np.savetxt('trumck_snrlumrad_cc_vst.txt',v_st)
    
            quan,unit = textstring(x,elem)
            grid[x].plot(tim,lum,colors[ind],lw=1.5,label=quan+'= %.1f '%(elem)+unit)
    # ax.axvline(t_st0,linestyle='--',lw=1.5)
            grid[x].set_xscale('log')
            grid[x].set_yscale('log')
    #ax.set_title('Mej = {}, n = {}, e51 = {}'.format(mej,nh,e51))
            grid[x].set_xlim(400,8.0e5)
            grid[x].set_ylim(1.0e20,1.0e24)
            grid[x].legend(fontsize='9')
            for axis in ['top','bottom','left','right']:
                grid[x].spines[axis].set_linewidth(1.2)
           
            
            grid[x+3].plot(tim,rad,colors[ind],lw=1.5,label=quan+'= %.1f '%(elem)+unit)
            grid[x+3].set_xscale('log')
            grid[x+3].set_yscale('log')
            grid[x+3].set_xlim(40,8.0e5)
            grid[x+3].set_ylim(0.3,80)
            grid[x+3].legend(loc=4,fontsize='9')
            for axis in ['top','bottom','left','right']:
                grid[x+3].spines[axis].set_linewidth(1.2)
            
            
            grid[x+6].plot(t_ed,v_ed,colors[ind])
            grid[x+6].plot(10**t_st,v_st,colors[ind],lw=1.5,label=quan+'= %.1f '%(elem)+unit)
            grid[x+6].set_xscale('log')
            grid[x+6].set_yscale('log')
            grid[x+6].set_xlim(400,8.0e5)
            grid[x+6].set_ylim(20,20000)
            grid[x+6].legend(loc=1,fontsize='9')
            for axis in ['top','bottom','left','right']:
                grid[x+6].spines[axis].set_linewidth(1.2)  
            

#py.tightlayout()
grid[0].set_ylabel('Radio Luminosity [ergs/s/Hz]',fontsize='13')
grid[3].set_ylabel('Radius [pc]',fontsize='13')
grid[6].set_ylabel('Ejecta Velocity [km/s]',fontsize='13')
grid[7].set_xlabel('Time [years]',fontsize='27')
py.savefig('indsnrgridplot.png')            
py.show()
