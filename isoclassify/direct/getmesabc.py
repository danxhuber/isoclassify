# interpolates bolometric correction from MIST tables, available at:
# http://waps.cfa.harvard.edu/MIST/model_grids.html

# available keys for photometric bands:
# Johnson BVRI: bc_b,bc_v,bc_r,bc_i
# Tycho BtVt: bc_bt,bc_vt
# Sloan ugriz:bc_us,bc_gs,bc_rs,bc_is,bc_zs
# 2MASS JHK: bc_j,bc_h,bc_k
# Kepler d51: bc_d51
# Gaia G: bc_ga'

# example:
# from getmesabc import *
# In [18]: from getmesabc import *
# In [19]: getbc(5777.,4.4,0.0,0.0,'bc_gs')
# Out[19]: -0.34045390880000004

import os
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import h5py
from isoclassify import DATADIR

def getbc(teff,logg,feh,av,band):

    fn = os.path.join(DATADIR, 'bcgrid.h5')
    bcmodel = h5py.File(fn, 'r')

    teffgrid = np.array(bcmodel['teffgrid'])
    logggrid = np.array(bcmodel['logggrid'])
    fehgrid = np.array(bcmodel['fehgrid'])
    avgrid = np.array(bcmodel['avgrid'])
    bcgrid = np.array(bcmodel[band])

    points = (teffgrid, logggrid, fehgrid, avgrid)
    values = bcgrid
    interp = RegularGridInterpolator(points, values) 

    xi = np.array([teff,logg,feh,av])
    bc = interp(xi)[0]
    return bc
