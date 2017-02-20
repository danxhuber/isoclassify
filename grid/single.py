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
    
    x.addspec([5065.,-99.0,-0.1],[120.,0.0,0.2])
    x.addseismo([231.,16.5],[10.,0.5])
    
    x.addjhk([6.025,5.578,5.496],[0.019,0.038,0.018])
    x.addplx(8.9536/1e3,0.7/1e3)
    
    paras=classify(input=x,model=model,dustmodel=0.,doplot=1)
