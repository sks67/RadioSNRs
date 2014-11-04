#Visualizing differential and cumulative histograms 

import numpy as np
import matplotlib.pyplot as plt
import os

def diam_cumhist(diam,f,d,R,fignum):
    path = '../Inputs/lmc_angdiams.txt'
    dist = 0.05 #Mpc
    lmcdiams_ang = np.loadtxt(path)
    lmcdiams_temp = ((lmcdiams_ang*np.pi)/(180.0*3600.0))*dist*(1.0e6)
    n1, bin1, patches1 = plt.hist(lmcdiams_temp,bins=10000,color='g',lw=3.0,histtype='step',cumulative=1)
    for ind,diams in enumerate(diam):
        diams = diams[np.nonzero(diams)]
        n2, bins2, patches2 = plt.hist(diams,bins=10000,histtype='step',alpha=0.2,color='purple',cumulative=1)
    plt.xlabel('Diameter[pc]',fontsize=15)
    plt.ylabel('N',fontsize=15)
    plt.title(r'$f=%.f,\ d=%.f pc,\ R=%.2e$'%(f,d,R),fontsize=15)
    ax = plt.gca()
    ax.tick_params(labelsize=15)
    plt.ylim(0,80)
    userdoc = os.path.join(os.getcwd(),'Plots')
    plt.xlim(8,1.7e2)
    plt.tight_layout()
    fig = plt.gcf()
    plt.savefig(os.path.join(userdoc,'DiamHist{}.png'.format(str(fignum))))
    plt.show()

