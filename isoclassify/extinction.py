import mwdust
import numpy as np
import pandas as pd
import astropy.units as units
from astropy.coordinates import SkyCoord
from dustmaps.bayestar import BayestarWebQuery

def query_dustmodel_coords(ra,dec):
    reddenMap = BayestarWebQuery(version='bayestar2017')
    sightLines = SkyCoord(ra*units.deg,dec*units.deg,frame='icrs')
    distanceSamples = np.loadtxt("isoclassify/data/distance-samples-green18.txt",delimiter=',')*1000. # In pc, from bayestar2017 map distance samples
    reddenContainer = reddenMap(sightLines,mode='best')
    
    del reddenMap # To clear reddenMap from memory
    
    dustModelDF = pd.DataFrame({'ra': [ra], 'dec': [dec]})

    for index in xrange(len(reddenContainer)):
        dustModelDF['av_'+str(round(distanceSamples[index],6))] = reddenContainer[index]
    
    return dustModelDF

def query_dustmodel_coords_allsky(ra,dec):
    reddenMap = mwdust.Combined15()
    sightLines = SkyCoord(ra*units.deg,dec*units.deg,frame='galactic')
    distanceSamples = np.loadtxt("isoclassify/data/distance-samples-green18.txt",delimiter=',')*1000. # Here this is just a container to keep the format of the code the same regardless of the utilized map. The combined map is *NOT* 3D, so the reddening should be the same at every distance.
    reddenContainer=reddenMap(sightLines.l.value,sightLines.b.value,distanceSamples/1000.)
    
    del reddenMap # To clear reddenMap from memory
    
    dustModelDF = pd.DataFrame({'ra': [ra], 'dec': [dec]})
    
    for index in xrange(len(reddenContainer)):
        dustModelDF['av_'+str(round(distanceSamples[index],6))] = reddenContainer[index]

    return dustModelDF

# R_lambda values to convert E(B-V) given by dustmaps to extinction in
# a given passband.  The two main caveats with this are: - strictly
# speaking only cardelli is consistent with the BC tables used in the
# MIST grid, but using wrong R_lambda's for the newer Green et
# al. dustmaps is (probably) worse.  - some values were interpolated
# to passbands that aren't included in the Schlafly/Green tables.

def extinction(law):
    if (law == 'cardelli'):
        out = {}
        with open("isoclassify/data/extinction-vector-cardelli.txt") as f:
            for line in f:
                (key,val) = line.split(',')
                out[key] = float(val)

    if (law == 'schlafly11'):
        out = {}
        with open("isoclassify/data/extinction-vector-schlafly11.txt") as f:
            for line in f:
                (key,val) = line.split(',')
                out[key] = float(val)

    if (law == 'schlafly16'):
        # see http://argonaut.skymaps.info/usage under "Gray Component". this is a lower limit.
        grayoffset=0.063
        out = {}
        with open("isoclassify/data/extinction-vector-schlafly16.txt") as f:
            for line in f:
                (key,val) = line.split(',')
                out[key] = float(val)+grayoffset
    return out
