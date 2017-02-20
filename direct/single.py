import sys
import ebf
import numpy as np
import scipy.interpolate
import pdb
import asfgrid
from astropy.io import ascii
from classify_direct import *


dnumodel = asfgrid.Seism()  
bcmodel = h5py.File('bcgrid.h5', 'r')
dustmodel = mwdust.Combined15()


x=obsdata()

x.addspec([5065.,-99.0,-0.1],[120.,0.0,0.2])
x.addseismo([231.,16.5],[10.,0.5])
paras=stparas(input=x,dnumodel=dnumodel,bcmodel=bcmodel,dustmodel=dustmodel,\
                      useav=0,dnucor=0,plot=0)

x.addjhk([6.025,5.578,5.496],[0.019,0.038,0.018])
x.addplx(8.9536/1e3,0.7/1e3)
paras=stparas(input=x,dnumodel=dnumodel,bcmodel=bcmodel,dustmodel=dustmodel,\
                      useav=0,dnucor=0,plot=0)

