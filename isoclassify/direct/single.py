import sys
import ebf
import numpy as np
import scipy.interpolate
import pdb
import asfgrid
import h5py
import mwdust

from astropy.io import ascii
from classify_direct import *



#dnumodel = asfgrid.Seism()  
bcmodel = h5py.File('bcgrid.h5', 'r')
dustmodel = mwdust.Combined15()


x=obsdata()

# example: HIP84970
x.addcoords(260.50241395,-24.99954638)
x.addmag([3.26],[0.01]) # magnitude for distance modulus; needs to be consistent with band keyword
x.addbv([3.03,3.26],[0.01,0.01])
x.addplx(7.48/1e3,0.17/1e3)

paras=stparas(input=x,bcmodel=bcmodel,dustmodel=dustmodel,\
                      useav=-99.,plot=1,band='v')


#raw_input(':')
#x.addspec([5065.,-99.0,-0.1],[120.,0.0,0.2])
#x.addseismo([231.,16.5],[10.,0.5])
