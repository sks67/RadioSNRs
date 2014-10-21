import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp

def interpvals(x,y,lik):
    X,Y = np.meshgrid(x,y)
    xnew = np.linspace(x.min(),x.max(),500)
    ynew = np.linspace(y.min(),y.max(),500)
    f = interp.interp2d(x,y,lik,kind='linear')
    znew = f(xnew,ynew)
    return (xnew,ynew,znew)

