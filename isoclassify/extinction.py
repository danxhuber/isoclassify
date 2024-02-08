import mwdust
import numpy as np
import pandas as pd
import astropy.units as units
from astropy.coordinates import SkyCoord
#from dustmaps.bayestar import BayestarWebQuery
import mwdust
from isoclassify import DATADIR

from . import PACKAGEDIR

def query_dustmodel_coords(ra,dec,dust):
    if dust == 'allsky':
        reddenMap = mwdust.Combined19()
        ext = extinction('green19')
    if dust == 'green19':
        reddenMap = mwdust.Green19()
        ext = extinction('green19')
    if dust == 'zero':
        reddenMap = mwdust.Zero()
        ext = extinction('cardelli')
    if dust == 'none':
        reddenMap = 0
        ext = extinction('cardelli')
        print('not dustmap, fitting for reddening if multiple photometric bands are provided.')
        return reddenMap,ext

    sightLines = SkyCoord(ra*units.deg,dec*units.deg,frame='icrs')
    sightLines = sightLines.transform_to('galactic')

    distanceSamples = np.loadtxt(f"{PACKAGEDIR}/data/distance-samples-green19.txt",delimiter=',')*1000.

    reddenContainer = reddenMap(sightLines.l.value,sightLines.b.value,distanceSamples/1000.)

    del reddenMap # To clear reddenMap from memory

    dustModelDF = pd.DataFrame({'ra': [ra], 'dec': [dec]})

    for index in range(len(reddenContainer)):
        dustModelDF['av_'+str(round(distanceSamples[index],6))] = reddenContainer[index]

    return dustModelDF,ext

# R_lambda values to convert E(B-V) given by dustmaps to extinction in
# a given passband.  The two main caveats with this are: - strictly
# speaking only cardelli is consistent with the BC tables used in the
# MIST grid, but using wrong R_lambda's for the newer Green et
# al. dustmaps is (probably) worse.  - some values were interpolated
# to passbands that aren't included in the Schlafly/Green tables.

def extinction(law):

    if (law == 'cardelli'):
        out = {}

        with open(f"{PACKAGEDIR}/data/extinction-vector-cardelli-iso.txt") as f:

            for line in f:
                (key,val) = line.split(',')
                out[key] = float(val)

    if (law == 'schlafly11'):
        out = {}

        with open(f"{PACKAGEDIR}/data/extinction-vector-schlafly11-iso.txt") as f:

            for line in f:
                (key,val) = line.split(',')
                out[key] = float(val)

    if (law == 'schlafly16'):
        # see http://argonaut.skymaps.info/usage under "Gray Component". this is a lower limit.
        grayoffset=0.063
        out = {}

        with open(f"{PACKAGEDIR}/data/extinction-vector-schlafly16-iso.txt") as f:

            for line in f:
                (key,val) = line.split(',')
                out[key] = float(val)+grayoffset

    if (law == 'green19'):
        out = {}

        with open(f"{PACKAGEDIR}/data/extinction-vector-green19-iso.txt") as f:

            for line in f:
                (key,val) = line.split(',')
                out[key] = float(val)
    return out
