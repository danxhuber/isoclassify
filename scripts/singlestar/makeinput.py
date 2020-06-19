import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import pdb
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astroquery.gaia import Gaia      
import astropy.units as units
from astropy.coordinates import SkyCoord
import astropy.units as u
import sys

#Vizier.VIZIER_SERVER = 'vizier.cfa.harvard.edu'

#gaiacols=Vizier(columns=['chi2AL','NgAL','Gmag','BP-RP','Plx','e_Plx','RPmag','BPmag','e_RPmag','e_BPmag'])
idstr=sys.argv[1]
feh=np.float(sys.argv[2])
fehe=np.float(sys.argv[3])
band=sys.argv[4]
dust=sys.argv[5]


simbad = Simbad.query_object(idstr)
coord = SkyCoord(ra=simbad['RA'], dec=simbad['DEC'], unit=(u.hourangle, u.deg))

plx=-99.
plxe=-99.
bpmags=-99.
bpmagse=-99.
rpmags=-99.
rpmagse=-99.
vmag=-99.
bmag=-99.
vmage=-99.
bmage=-99.
jmags=-99.
hmags=-99.
kmags=-99.
jmagse=-99.
hmagse=-99.
kmagse=-99.
ra=-99.
dec=-99.

print('querying Gaia ...')
try:
    width = u.Quantity(0.01, u.deg)     
    height = u.Quantity(0.01, u.deg)     
    gaiaq = Gaia.query_object_async(coordinate=coord[0], width=width, height=height) 
    plx=gaiaq['parallax'][0]+0.05
    plxe=gaiaq['parallax_error'][0]
    bpmags=gaiaq['phot_bp_mean_mag'][0]
    bpmagse=0.01
    rpmags=gaiaq['phot_rp_mean_mag'][0]
    rpmagse=0.01
    print('done')
except:
    gaiaq = 0
    print('failed')

print('querying Tycho ...')
tychocols=Vizier(columns=['RAmdeg','DEmdeg','VTmag','BTmag','e_VTmag','e_BTmag'])
try:
    tycho=tychocols.query_object(idstr,catalog="I/259",radius=5.0*units.arcsec) 
    vmag=tycho[0]['VTmag'][0]
    bmag=tycho[0]['BTmag'][0]
    vmage=tycho[0]['e_VTmag'][0]
    bmage=tycho[0]['e_BTmag'][0]
    print('done')
except:
    tycho = 0
    print('failed')

print('querying 2MASS ...')
try:
    tmass = Vizier.query_object(idstr,catalog="II/246",radius=5.0*units.arcsec)
    jmags=tmass[0]['Jmag'][0]
    hmags=tmass[0]['Hmag'][0]
    kmags=tmass[0]['Kmag'][0]
    jmagse=tmass[0]['e_Jmag'][0]
    hmagse=tmass[0]['e_Hmag'][0]
    kmagse=tmass[0]['e_Kmag'][0]
    ra=tmass[0]['RAJ2000'][0]
    dec=tmass[0]['DEJ2000'][0] 
    print('done')
except:
    tmass=0
    print('failed')
    
ascii.write([[idstr],[ra],[dec],[vmag],[vmage],[bmag],[bmage],[bpmags],[bpmagse],[rpmags],[rpmagse],[jmags],[jmagse],[hmags],[hmagse],[kmags],[kmagse],[plx/1e3],[plxe/1e3],[band],[dust],[-99.],[-99.],[feh],[fehe]],'input_direct.csv',names=['id_starname','ra','dec','vtmag','vtmag_err','btmag','btmag_err','bpmag','bpmag_err','rpmag','rpmag_err','jmag','jmag_err','hmag','hmag_err','kmag','kmag_err','parallax','parallax_err','band','dust','logg','logg_err','feh','feh_err'],delimiter=',',overwrite=True)
