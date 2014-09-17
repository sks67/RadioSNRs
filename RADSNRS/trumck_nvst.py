#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Same as trumck_checksnaps_easier.py, except here I attempt 
# several optimization methods.
#
# NOTE: If program works, delete trumck_checksnaps_easier.py 
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
import os
import datetime as dt
import numpy as np
import matplotlib.pyplot as py
import math as mt
import scipy.interpolate as interp
import scipy.integrate as integ
from astropy.io import ascii
from scipy import stats
import ConfigParser

Config = ConfigParser.ConfigParser()
Config.read('config.ini')
nlmc = Config.getint('ObsParameters','numberofsnrs')
dist = Config.getfloat('ObsParameters','distance')   #Mpc
path_to_file =  Config.get('InputFiles','path') 
lumfile =  Config.get('InputFiles','luminosity') 

lmclums_jy = np.loadtxt(path_to_file+lumfile)
lmclums=(1.0e24)*(lmclums_jy*1000.)*1.2*dist*dist


pc = 3.086e18



#---------------------------------------------------------------#                                                                                                                  
#------------PICKING RANDOM NH IN A COLUMN DENSITY--------------#                                                                                                                   
#---------------------------------------------------------------#                                                                                                                  
def poissonProb(j,n):
    return (np.exp(-n)*(n**j))/mt.factorial(j)

def mostprobpval(kvals):
    numbins=5
    probs,bin_edges = np.histogram(np.log10(kvals),bins=numbins,density=True)
    arg= np.argmax(probs)
    return 10**(bin_edges[arg]+((bin_edges[arg+1]-bin_edges[arg])/2.))

def randomnh(rho_col,z0):
    #Defining the PDF
    mean=0.0   #pc
    stdev=z0/np.sqrt(2)
    z = np.random.normal(loc=mean,scale=stdev,size=1)
    n_z = (rho_col/(pc*z0*np.sqrt(3.14)))*np.exp(-(z**2)/(z0**2))
    return n_z

def radiolightcurve(lmcthick,nh,tborn,ejmas,energ,nprof):
    #thick_lim = 5*(lmcthick/pc)*0.5
    tsnap_array = np.logspace(3.0,7.0,500)
    histlum_array = np.zeros_like(tsnap_array)
    nsnrs = np.zeros(tsnap_array.size)
    
#Declaring constants for Luminosity Calculation
    epsb = 0.01
    alpha = 1.0 
    epse = epsb*alpha
    pp = 3.0
    compf = 4.0
    ff = 0.38
    nei = 1.14
    cl = 6.27e18
    c5 = 7.52e-24
    c6 = 7.97e-41
    dist = 50*1000*3.086e18
    nu = 1.4e9
    

    for xx in range(tsnap_array.size):
        tsnap = tsnap_array[xx] #years
#---------------------------------------------------------------#                                                                                                                  
#---------------------------------------------------------------#                                                                                                                  
#---------------------------------------------------------------#   

#random locations of densities and times for snrate = 2.0e3
               
        lum = np.zeros(nh.size) #800 is the size of lluni2 for the sedov taylor phases
        tim = np.zeros(tborn.size) #800 is the size of lluni2 for the sedov taylor phases
        for i in range(nh.size):
            if tborn[i]>tsnap:
                break
            nt2=800
            lluni2 = 0.0
            mu=1.4
            n0 = randomnh(nh[i],lmcthick) #dimensionless, multiplied by 1cm^3
            mp = 1.67e-24
            ismrho = n0*mu*mp #in units of cm^-3
            e51 = energ[i] #in units of 10^51 ergs of energy released per sne
            mej = ejmas[i] #in units of solar masses
            if n0==0.0:
                print n0, nh[i]
#Characteristic scales
            tch = 423*e51*(mej**(5.0/6.0))*(n0**(-1.0/3.0)) #years
            rch = 3.07*(mej**(1.0/3.0))*(n0**(-1.0/3.0)) #pcs
            vch = 7090*(e51**0.5)*(mej**(-0.5)) #km/s
                
#variables 
   
            t_ed = tsnap-tborn[i]
        #if t_ed<0.0:
         #   continue
            if(nprof[i]==7.0):
                tstar_st = 0.732
                t_st0 = tstar_st*tch
                tstar_ed = t_ed/tch
                vstar_ed = 0.606*tstar_ed**(-(3.0/7.0)) if t_ed<t_st0 else 0.569*((1.42*tstar_ed - 0.312)**(-3.0/5.0))
                rstar_ed = 1.06*tstar_ed**(4.0/7.0) if t_ed<t_st0 else (1.42*tstar_ed - 0.312)**(2.0/5.0)
                v_ed = vstar_ed*vch  #in km/s                                                                                                                                         
                r_ed = rstar_ed*rch  #in par
            
#--------------------------------------------------------#    
            elif(nprof[i]==12.0):
                tstar_st = 0.424
                t_st0 = tstar_st*tch
                tstar_ed = t_ed/tch
                vstar_ed = 0.545*tstar_ed**(-(3.0/7.0)) if t_ed<t_st0 else 0.569*((1.42*tstar_ed - 0.28)**(-3.0/5.0))
                rstar_ed = 0.953*tstar_ed**(4.0/7.0) if t_ed<t_st0 else (1.42*tstar_ed - 0.28)**(2.0/5.0)
                v_ed = vstar_ed*vch  #in km/s                                                                                                                                         
                r_ed = rstar_ed*rch  #in par

            vrad = 200.
            if (v_ed<=vrad):
                lum[i]=0.0
                tim[i]=t_ed
                continue
#----------------------------------------------------------------------------------------------#                                                                                    #-------------------------------ED STAGE LUMINOSITY AND RADIUS---------------------------------#                        
#----------------------------------------------------------------------------------------------#                                                                                    
            rho0 = ismrho/((1.0e-24)*mu*1.67)
            bism = 9e-6*rho0**0.47    #magnetic field for galaxies scales as 0.47 (Krutcher et al. 1999)                                                                              
            rr1 = r_ed*(3.086e18) #LAURA - cm                                                                                                                                          
            vv1 = v_ed*1.0e5 #LAURA - cm/s                                                                                                                                              
            bb0 = np.sqrt(8.0*mt.pi*epsb*ismrho*vv1*vv1)
            bb1 = np.sqrt(8.0*mt.pi*epsb*ismrho*vv1*vv1)
            if bb1<(4.0*bism):
                bb1=4.0*bism
            me = 9.11e-28  #in grams
            c = 3.0e10
            gammam = (mu*epse*(pp-2)*((vv1/c)**2)*mp)/((pp-1)*compf*nei*me)
            if gammam<1.0:
                gammam=1.0
   
            elow = gammam*me*c**2  #Eq. 10 Chevalier, 98                                                                                                                       
            n0 = (alpha*(bb0**2)*(pp-2)*(elow**(pp-2.0)))/(8.0*mt.pi)
            ss11 = (4.0/3.0)*ff*rr1
            nu11 = 2.0*cl*((ss11*c6*n0)**(2.0/(pp+4.0)))*(bb1**((pp+2.0)/(pp+4.0)))
            ss1 = (c5/c6)*(bb1**(-0.5))*((nu11/(2.0*cl))**2.5)
            jj1 = ((nu/nu11)**2.5)*(1.0-np.exp(-1.0*(nu/nu11)**(-1.0*(pp+4.0)/2.0)))
            f_nu1 = (ss1*jj1*mt.pi*rr1**2)/(1.0e-26*dist**2)
            lluni = f_nu1*1.0e-26*4.0*mt.pi*dist**2
            if (lluni<=9.0e21):
                lluni=0.0
        
            lum[i]=lluni
            tim[i]=t_ed

        histlum = lum[np.nonzero(lum)]
        nsnrs[xx]=histlum.size
        if xx%10==0:
            print tsnap,'\t',nsnrs[xx] 
#-------------------------------------------------------------------#                                                                                     

      
   
#-------------------------------------------------------------------#
#-------------- LIKELIHOOD TESTING (Badenes2010) -------------------#
#-------------------------------------------------------------------#


#Saving stuff to the IpythonStuff folder for analysis    
    userdoc = os.path.join(os.getcwd(),'DataAnalysis')    
    #np.savetxt(os.path.join(userdoc,'paramtst_lumhist.txt'),histlum_array)
    np.savetxt(os.path.join(userdoc,'paramtst_nsnrs.txt'),nsnrs)
    
    choice=0.0
    if choice==1.0:
#CALCULATE N(OBS) AND N(MODEL) PER BIN
    
        l_cutoff = 2.0e23 #ergs/s/Hz
        
        obs_bins = np.logspace(np.log10(l_cutoff),np.log10(lmclums.max()),6)
        n,bins =np.histogram(lmclums,bins=obs_bins)
        n2 = np.zeros((histlum_array.shape[0],obs_bins.size-1))
        for ind,lums in enumerate(histlum_array):
            lums = lums[np.nonzero(lums)]
            n2[ind], bins2 = np.histogram(lums,bins=obs_bins)
        
#CALCULATING LIKELIHOOD
        avg_n = np.array([np.mean(col) for col in np.transpose(n2)])
        likhood_temp = np.array([poissonProb(n[ind],avg_n[ind]) for ind in range(n.size)])
        likhood = np.prod(likhood_temp)
        return likhood
    

