import pdb
import numpy as np
import matplotlib.pyplot as plt
import fnmatch
from scipy.ndimage.filters import gaussian_filter
import time
from scipy.interpolate import interp1d
import pandas as pd

def binpdf(x,y,step,iname,dustmodel):
    if np.max(x) - np.min(x) < step: # Ensure step sizes are smaller than the input array value differences
        step = (np.max(x) - np.min(x))/10.
        print('Default '+iname+' step size too large. Using ' + str(step) + ' step size. Treat results with caution.')

    xax = np.arange(np.min(x),np.max(x),step)

    if fnmatch.fnmatch(iname,'*age*'):
        xax = 10.**(np.arange(71.)/(70./2.15) - 1.) # Exact grid ages (log + linear)
        xax = np.concatenate((xax[:50],xax[50]+np.arange(67)*0.25))

    elif ( (isinstance(dustmodel,pd.DataFrame) == False) & (fnmatch.fnmatch(iname,'*avs*'))):
        grid=np.unique(x)
        spacing=grid[1]-grid[0]
        xax = np.arange(grid[0]-spacing/4.,grid[len(grid)-1]+spacing/4.,spacing)
        step=spacing

    else:
        xax= xax+step/2.

    yax = np.zeros(len(xax))

    digitized = np.digitize(x, xax)
    yax = [y[digitized == i].sum() for i in range(1, len(xax)+1)]

    '''
    for r in range(0,len(xax)-1):
    ix = np.where((x > xax[r]) & (x <= xax[r+1]))
    yax[r] = np.sum(y[ix[0]])

    plt.clf()

    plt.plot(xax,yax2/np.max(yax2))
    plt.plot(xax,yax/np.max(yax))
    pdb.set_trace()
    '''
    #yax = gaussian_filter(yax,1.5)
    #if fnmatch.fnmatch(iname,'*avs*'):
    #         pdb.set_trace()
    #pdb.set_trace()
    yax = yax/np.sum(yax)

    return xax,yax


def getstat(xax,yax):
    cdf = np.cumsum(yax)
    # bad hack, needs to be better
    if (np.min(cdf) > 0.16):
        cdf[np.argmin(cdf)]=0.16
    if (np.max(cdf) < 0.84):
        cdf[np.argmax(cdf)]=0.84
    #pdb.set_trace()
    ppf = interp1d(cdf,xax) # percent point function
    p16, med, p84 = ppf([0.16,0.50,0.84])
    emed1  = med - p16
    emed2  = p84 - med
    return med,emed2,emed1


def getpdf(x,y,step,fixed,name,dustmodel):
    if fixed == 0:
        pos=np.argmax(y)
        steps=x[pos]*step
    else:
        steps=step

    xax,yax = binpdf(x=x,y=y,step=steps,iname=name,dustmodel=dustmodel)
    med,emed1,emed2 = getstat(xax,yax)

    if ( (isinstance(dustmodel,pd.DataFrame) == False) & (fnmatch.fnmatch(name,'*avs*'))):
        return xax,yax,med,emed1,emed2
    if fnmatch.fnmatch(name,'*feh*'):
        return xax,yax,med,emed1,emed2
    if fnmatch.fnmatch(name,'*age*'):
        return xax,yax,med,emed1,emed2

    newstep = ((emed1+emed2)/2.)/10.
    if newstep > 0.:
        steps=newstep

    xax,yax = binpdf(x=x,y=y,step=steps,iname=name,dustmodel=dustmodel)
    med,emed1,emed2 = getstat(xax,yax)

    if fnmatch.fnmatch(name,'*rho*'):
        xax=10**xax
        #yax=10**yax
        med,emed1,emed2 = getstat(xax,yax)

    if fnmatch.fnmatch(name,'*lum*'):
        xax=10**xax
        #yax=10**yax
        med,emed1,emed2 = getstat(xax,yax)

    #if plot > 0:
    #        plt.subplot(7,2,plot)
    #        plt.plot(xax,np.cumsum(yax))
    #        plt.subplot(7,2,plot+1)
    #        plt.plot(xax,yax)

    return xax,yax,med,emed1,emed2
