#%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
from classify_grid import *
import os, ebf
from astropy.io import ascii
import time
#import mwdust

if __name__ == '__main__':

    homedir=os.path.expanduser('~/')
    model=ebf.read(homedir+'science/models/MIST/mesa.ebf')

    model['rho']=np.log10(model['rho'])
    # do this to turn off scaling relation corrections
    model['fdnu'][:]=1.
    model['avs']=np.zeros(len(model['teff']))
    model['dis']=np.zeros(len(model['teff']))

    #x.addcoords(338.3683920,-9.0227690)
    #dustmodel = mwdust.Combined15()

    x=obsdata()
    
    x.addspec([5801.,-99.0,-0.07],[80.,0.0,0.1])
    x.addseismo([1240.,63.5],[70.,1.5])
    x.addjhk([10.369,10.07,10.025],[0.022,0.018,0.019])
    x.addgriz([11.776,11.354,11.238,11.178],[0.02,0.02,0.02,0.02])

    #x.addplx(2.71/1e3,0.08/1e3)
    paras=classify(input=x,model=model,dustmodel=0.,doplot=1)