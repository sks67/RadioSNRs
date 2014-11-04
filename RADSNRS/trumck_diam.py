
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  NOTE: IF PROGRAM WORKS, DELETE trumck_diam.py
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
diamfile = Config.get('InputFiles','size')
diam = np.loadtxt(path_to_file+diamfile)  #arcsecs                                                                                                     
lmcdiams = ((diam*mt.pi)/(180.0*3600.0))*dist*(1.0e6) #pc


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
    tsnap_array = np.linspace(1.8e6,2.4e6,50)
    diam_array = np.zeros((tsnap_array.size,9000))
    nsnrs = np.zeros(tsnap_array.size)
    
#Declaring constants for Luminosity Calculation
    epsb = 0.01
    alpha = 1.0 
    epse = epsb*alpha
    pp = 2.5
    compf = 4.0
    ff = 0.38
    nei = 1.14
    cl = 6.27e18
    c5 = 9.68e-24
    c6 = 8.10e-41
    dist = 50*1000*3.086e18
    nu = 1.4e9
    

    for xx in range(tsnap_array.size):
        tsnap = tsnap_array[xx] #years
#---------------------------------------------------------------#                                                                                     
#---------------------------------------------------------------#                                                                                       
#---------------------------------------------------------------#   

#random locations of densities and times for snrate = 2.0e3
               
        rad = np.zeros(nh.size)
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
            tch = 423*(e51**(-0.5))*(mej**(5.0/6.0))*(n0**(-1.0/3.0)) #years
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
                tim[i]=t_ed
                rad[i]=0.0
                continue
        
            tim[i]=t_ed
            rad[i]=r_ed

        histrad = rad[np.nonzero(rad)]
        diam_array[xx,0:histrad.size]=2.0*histrad
#-------------------------------------------------------------------#                                                                                     

      
   
#-------------------------------------------------------------------#
#-------------- LIKELIHOOD TESTING (Badenes2010) -------------------#
#-------------------------------------------------------------------#


#Saving stuff to the IpythonStuff folder for analysis    
    userdoc = os.path.join(os.getcwd(),'DataAnalysis')    
    np.savetxt(os.path.join(userdoc,'paramtst_diamhist.txt'),diam_array)
 
#CALCULATE N(OBS) AND N(MODEL) PER BIN
    
    diam_cutoff = 80.0 #pc

    obs_bins = np.linspace(0,diam_cutoff,6)
    n,bins =np.histogram(lmcdiams,bins=obs_bins)
    n2 = np.zeros((diam_array.shape[0],obs_bins.size-1))
    for ind,lums in enumerate(diam_array):
        lums = lums[np.nonzero(lums)]
        n2[ind], bins2 = np.histogram(lums,bins=obs_bins)
    

        
#CALCULATING LIKELIHOOD
    avg_n = np.mean(n2,axis=0)
    likhood_temp = np.array([poissonProb(n[ind],avg_n[ind]) for ind in range(n.size)])
    likhood = np.prod(likhood_temp)
    return likhood
    

