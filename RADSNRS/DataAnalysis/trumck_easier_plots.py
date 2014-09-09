import numpy as np
import matplotlib.pyplot as py
import subprocess
import scipy.interpolate as interp
import scipy.integrate as integ
from mpl_toolkits.axes_grid1 import Grid

inp = np.genfromtxt('trumck_inputs.txt',dtype='str',comments='#',skip_header=4)

histlum_array = np.loadtxt('histlum_array.txt')
histype_array = np.loadtxt('histype_array.txt')
histrho_array = np.loadtxt('histrho_array.txt')
kstest_array = np.loadtxt('kstest_array')
size_array = np.loadtxt('size_array')
lmcrho = np.loadtxt(inp[-1])
lmclums = np.loadtxt('paramtst_lumlmc.txt')

frac_array = np.array([0.2,0.4,0.6,0.8,1.0])
thick_array = np.array([10,20,30,40,50])

fig = py.figure(1,figsize=(12,12))
grid = Grid(fig,rect=111,nrows_ncols=(thick_array.size,frac_array.size),axes_pad=0.0,label_mode='L')
fig2 = py.figure(2,figsize=(10,10))
grid2 = Grid(fig2,rect=111,nrows_ncols=(thick_array.size,frac_array.size),axes_pad=0.0,label_mode='L')

x=0
for j in range(thick_array.size):
    for i in range(frac_array.size):
         #Plot lum-rho grid
        print x
        histlum = histlum_array[x][np.nonzero(histlum_array[x])]
        histype = histype_array[x][np.nonzero(histype_array[x])]
        histrho = histrho_array[x][np.nonzero(histrho_array[x])]
        
        lum_1a = histlum[np.where(histype==7.0)[0]]
        lum_cc = histlum[np.where(histype==12.0)[0]]
        rho_1a = histrho[np.where(histype==7.0)[0]]
        rho_cc = histrho[np.where(histype==12.0)[0]]
        grid[x].plot(lmcrho,lmclums,'ko')
        grid[x].plot(rho_1a,lum_1a,'bo',mew=1.3,mec='b',mfc='none')
        grid[x].plot(rho_cc,lum_cc,'ro',mew=1.3,mec='r',mfc='none')
        grid[x].set_xscale('log')
        grid[x].set_yscale('log')
        grid[x].set_xlim(2.0e19,3.0e23)
        grid[x].set_ylim(5.0e22,2.0e26)
        textstr = 'p = %.2e\nN = %.f\n$\mathbf{z_0 =}$ %.f pc'%(kstest_array[x],size_array[x],thick_array[j])
       # props = dict(facecolor='None')
        grid[x].text(0.05,0.95,textstr,transform=grid[x].transAxes,fontsize=11,verticalalignment='top')


        n2,bins2,patches2=grid2[x].hist(histlum,10000,histtype='step',normed=False,cumulative=-1)
        n1,bins1,patches1=grid2[x].hist(lmclums,10000,histtype='step',normed=False,cumulative=-1)
        grid2[x].set_xscale('log')
        grid2[x].set_xlim(5.0e22,8.0e25)
        grid2[x].set_ylim(0,79)
        grid2[x].get_yaxis().set_ticks([10,30,50,70])
        grid2[x].axvline(x=3.0e23,color='k',linestyle='--',lw=1.5)
        patches2[0].set_xy(patches2[0].get_xy()[1:])
        patches1[0].set_xy(patches1[0].get_xy()[1:])
       # props = dict(boxstyle='square',facecolor='None')
        grid2[x].text(0.4,0.95,textstr,transform=grid2[x].transAxes,fontsize=11,verticalalignment='top')
        print frac_array[i],'\t',kstest_array[x]
        x=x+1


grid[22].set_xlabel(r'$\rho_{col}\ [cm^{-2}]$',fontsize=30)                                                                                                              
grid[10].set_ylabel(r'1.4 GHz Luminosity [ergs/s/Hz]',fontsize=30)
grid2[10].set_ylabel(r'N(>L)',fontsize=24)
grid2[22].set_xlabel(r'1.4 GHz Luminosity [ergs/s/Hz]',fontsize=24)

for i in range(5):
    textstr = 'f = %.2f'%(frac_array[i])
    grid[i].set_title(textstr,fontsize=17,weight='bold')
    grid2[i].set_title(textstr,fontsize=17,weight='bold')

np.savetxt('kstest_array',kstest_array)
np.savetxt('size_array',size_array)
fig.tight_layout()
fig2.tight_layout()
fig.savefig('gridplot_lumrad.png')
fig2.savefig('gridplot_lumhist.png')
py.close()        
