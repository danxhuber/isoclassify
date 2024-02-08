import copy
import time,pdb

import ephem
import pandas as pd
import numpy as np
from astropy.io import ascii
from itertools import product

from .pdf import *  # part of isoclassify package (to do make explicit import)
from .priors import * # part of isoclassify package (to do make explicit import)
from .plot import * # part of isoclassify package (to do make explicit import)

class obsdata():
    def __init__(self):
        self.plx = -99.0
        self.plxe = -99.0

        self.teff = -99.0
        self.teffe = -99.0
        self.logg = -99.0
        self.logge = -99.0
        self.feh = -99.0
        self.fehe = -99.0

        self.lum = -99.0
        self.lume = -99.0

        self.bmag = -99.0
        self.bmage = -99.0
        self.vmag = -99.0
        self.vmage = -99.0

        self.btmag = -99.0
        self.btmage = -99.0
        self.vtmag = -99.0
        self.vtmage = -99.0

        self.dmag = -99.0
        self.dmage = -99.0

        self.gmag = -99.0
        self.gmage = -99.0
        self.rmag = -99.0
        self.rmage = -99.0
        self.imag = -99.0
        self.image = -99.0
        self.zmag = -99.0
        self.zmage = -99.0
        self.jmag = -99.0
        self.jmage = -99.0
        self.hmag = -99.0
        self.hmage = -99.0
        self.kmag = -99.0
        self.kmage = -99.0

        self.gamag = -99.0
        self.gamage = -99.0
        self.bpmag = -99.0
        self.bpmage = -99.0
        self.rpmag = -99.0
        self.rpmage = -99.0

        self.numax = -99.0
        self.numaxe = -99.0
        self.dnu = -99.0
        self.dnue = -99.0

    def addspec(self,value,sigma):
        self.teff = value[0]
        self.teffe = sigma[0]
        self.logg = value[1]
        self.logge = sigma[1]
        self.feh = value[2]
        self.fehe = sigma[2]

    def addlum(self,value,sigma):
        self.lum = value[0]
        self.lume = sigma[0]

    def addbv(self,value,sigma):
        self.bmag = value[0]
        self.bmage = sigma[0]
        self.vmag = value[1]
        self.vmage = sigma[1]

    def addbvt(self,value,sigma):
        self.btmag = value[0]
        self.btmage = sigma[0]
        self.vtmag = value[1]
        self.vtmage = sigma[1]

    def addgriz(self,value,sigma):
        self.gmag = value[0]
        self.gmage = sigma[0]
        self.rmag = value[1]
        self.rmage = sigma[1]
        self.imag = value[2]
        self.image = sigma[2]
        self.zmag = value[3]
        self.zmage = sigma[3]

    def addjhk(self,value,sigma):
        self.jmag = value[0]
        self.jmage = sigma[0]
        self.hmag = value[1]
        self.hmage = sigma[1]
        self.kmag = value[2]
        self.kmage = sigma[2]

    def addgaia(self,value,sigma):
        self.gamag = value[0]
        self.gamage = sigma[0]
        self.bpmag = value[1]
        self.bpmage = sigma[1]
        self.rpmag = value[2]
        self.rpmage = sigma[2]

    def addplx(self,value,sigma):
        self.plx = value
        self.plxe = sigma

    def adddmag(self,value,sigma):
        self.dmag = value
        self.dmage = sigma

    def addseismo(self,value,sigma):
        self.numax = value[0]
        self.numaxe = sigma[0]
        self.dnu = value[1]
        self.dnue = sigma[1]

    def addcoords(self,value1,value2):
        self.ra = value1
        self.dec = value2

class resdata():
    def __init__(self):
        self.teff = 0.0
        self.teffep = 0.0
        self.teffem = 0.0
        self.teffpx = 0.0
        self.teffpy = 0.0
        self.logg = 0.0
        self.loggep = 0.0
        self.loggem = 0.0
        self.loggpx = 0.0
        self.loggpy = 0.0
        self.feh = 0.0
        self.fehep = 0.0
        self.fehem = 0.0
        self.fehpx = 0.0
        self.fehpy = 0.0
        self.rad = 0.0
        self.radep = 0.0
        self.radem = 0.0
        self.radpx = 0.0
        self.radpy = 0.0
        self.mass = 0.0
        self.massep = 0.0
        self.massem = 0.0
        self.masspx = 0.0
        self.masspy = 0.0
        self.rho = 0.0
        self.rhoep = 0.0
        self.rhoem = 0.0
        self.rhopx = 0.0
        self.rhopy = 0.0
        self.lum = 0.0
        self.lumep = 0.0
        self.lumem = 0.0
        self.lumpx = 0.0
        self.lumpy = 0.0
        self.age = 0.0
        self.ageep = 0.0
        self.ageem = 0.0
        self.agepx = 0.0
        self.agepy = 0.0
        self.avs = 0.0
        self.avsep = 0.0
        self.avsem = 0.0
        self.avspx = 0.0
        self.avspy = 0.0
        self.dis = 0.0
        self.disep = 0.0
        self.disem = 0.0
        self.dispx = 0.0
        self.dispy = 0.0

        self.teffsec = 0.0
        self.teffsecep = 0.0
        self.teffsecem = 0.0
        self.teffsecpx = 0.0
        self.teffsecpy = 0.0
        self.radsec = 0.0
        self.radsecep = 0.0
        self.radsecem = 0.0
        self.radsecpx = 0.0
        self.radsecpy = 0.0
        self.loggsec = 0.0
        self.loggsecep = 0.0
        self.loggsecem = 0.0
        self.loggsecpx = 0.0
        self.loggsecpy = 0.0
        self.rhosec = 0.0
        self.rhosecep = 0.0
        self.rhosecem = 0.0
        self.rhosecpx = 0.0
        self.rhosecpy = 0.0
        self.masssec = 0.0
        self.masssecep = 0.0
        self.masssecem = 0.0
        self.masssecpx = 0.0
        self.masssecpy = 0.0

class extinction():
    def __init__(self):
        self.ab = 1.3454449
        self.av = 1.00

        self.abt = 1.3986523
        self.avt = 1.0602271

        self.ag = 1.2348743
        self.ar = 0.88343449
        self.ai = 0.68095687
        self.az = 0.48308430

        self.aj = 0.28814896
        self.ah = 0.18152716
        self.ak = 0.11505195

        self.aga=1.2348743


def classify(input, model, dustmodel=0, plot=1, useav=-99.0, ext=-99.0, band=''):
    """
    Run grid based classifier

    Args:
        input (object): input object
        model (dict): dictionary of arrays
        dustmodel (Optional[DataFrame]): extinction model
        useav (float):
        ext (float):
    """

    ## constants
    gsun = 27420.010
    numaxsun = 3090.0
    dnusun = 135.1
    teffsun = 5772.0

    # bolometric correction error; kinda needs to be motivated better ...
    bcerr = 0.03

    ## extinction coefficients
    extfactors = ext

    ## class containing output results
    result = resdata()

    # calculate colors + errors:
    bvcol = input.bmag - input.vmag
    bvtcol = input.btmag - input.vtmag
    grcol = input.gmag - input.rmag
    ricol = input.rmag - input.imag
    izcol = input.imag - input.zmag
    gicol = input.gmag - input.imag
    rzcol = input.rmag - input.zmag
    gzcol = input.gmag - input.zmag
    jhcol = input.jmag - input.hmag
    hkcol = input.hmag - input.kmag
    jkcol = input.jmag - input.kmag
    bpgacol = input.bpmag - input.gamag
    garpcol = input.gamag - input.rpmag
    bprpcol = input.bpmag - input.rpmag
    vjcol = input.vmag - input.jmag
    vtjcol = input.vtmag - input.jmag
    gjcol = input.gmag - input.jmag
    rjcol = input.rmag - input.jmag
    vkcol = input.vmag - input.kmag
    vtkcol = input.vtmag - input.kmag
    gkcol = input.gmag - input.kmag
    rkcol = input.rmag - input.kmag
    gajcol = input.gamag - input.jmag
    gakcol = input.gamag - input.kmag

    bvcole = np.sqrt(input.bmage**2 + input.vmage**2)
    bvtcole = np.sqrt(input.btmage**2 + input.vtmage**2)
    grcole = np.sqrt(input.gmage**2 + input.rmage**2)
    ricole = np.sqrt(input.rmage**2 + input.image**2)
    izcole = np.sqrt(input.image**2 + input.zmage**2)
    gicole = np.sqrt(input.gmage**2 + input.image**2)
    rzcole = np.sqrt(input.rmage**2 + input.zmage**2)
    gzcole = np.sqrt(input.gmage**2 + input.zmage**2)
    jhcole = np.sqrt(input.jmage**2 + input.hmage**2)
    hkcole = np.sqrt(input.hmage**2 + input.kmage**2)
    jkcole = np.sqrt(input.jmage**2 + input.kmage**2)
    bpgacole = np.sqrt(input.bpmage**2 + input.gamage**2)
    garpcole = np.sqrt(input.gamage**2 + input.rpmage**2)
    bprpcole = np.sqrt(input.bpmage**2 + input.rpmage**2)
    vjcole = np.sqrt(input.vmage**2 + input.jmage**2)
    vtjcole = np.sqrt(input.vtmage**2 + input.jmage**2)
    gjcole = np.sqrt(input.gmage**2 + input.jmage**2)
    rjcole = np.sqrt(input.rmage**2 + input.jmage**2)
    vkcole = np.sqrt(input.vmage**2 + input.kmage**2)
    vtkcole = np.sqrt(input.vtmage**2 + input.kmage**2)
    gkcole = np.sqrt(input.gmage**2 + input.kmage**2)
    rkcole = np.sqrt(input.rmage**2 + input.kmage**2)
    gajcole = np.sqrt(input.gamage**2 + input.jmage**2)
    gakcole = np.sqrt(input.gamage**2 + input.kmage**2)

    # Compute extra color error term based on underestimation of stellar teff errors with nominal 2% error floor:
    if ((input.gmag > -99.0) & (input.kmag > -99.0)):
        gkexcole = compute_extra_gk_color_error(gkcol)
	# Determine which gK error term is greater and use that one:
        print("g - K error from photometry: ",gkcole)
        print("g - K error from best-fit polynomial: ",gkexcole)
        gkcole = max(gkcole,gkexcole)
        print("Using g - K error: ",gkcole)

    # apparent mag to use for distance estimation. set by "band" input
    redmap = -99.0 
    #if (getattr(input,band) > -99.):
    if pd.notnull(band):
        redmap = getattr(input,band)
        redmape = getattr(input,band+'e')
        model_mabs = model[band]
        # correct for companion
        if (input.dmag != -99.):
            dx=-0.4*input.dmag
            dxe=-0.4*input.dmage
            cor=2.5*np.log10(1.+10**dx)
            redmap = redmap+cor
            redmape = np.sqrt( redmape**2 + (dxe*2.5*10**dx/(1.+10**dx))**2)

    # absolute magnitude
    if (input.plx > -99.0):
        mabs = -5.0 * np.log10(1.0 / input.plx) + redmap + 5.0
        mabse = np.sqrt(
            (-5.0 / (input.plx * np.log(10)))**2 * input.plxe**2
            + redmape**2 + bcerr**2)

        # Also compute extra error term for M-dwarfs with K band mags only:
        if (mabs > 4.0) and (input.kmag > -99.0):
            print("M-dwarf with K band magnitude detected!")
            mabseex = compute_extra_MK_error(mabs)
            print("M_K from photometry: ",mabse)
            print("M_K error from best-fit polynomial: ",mabseex)
            mabse = np.sqrt(mabse**2 + mabseex**2)
            print("After adding in quadrature, using M_K error: ",mabse)
    else:
        mabs = -99.0
        mabse = -99.0


    # pre-select model grid; first only using reddening-independent quantities
    sig = 4.0
    um = np.arange(0,len(model['teff']),1)

    if (input.teff > -99.0):
        ut=np.where((model['teff'] > input.teff-sig*input.teffe) & \
        (model['teff'] < input.teff+sig*input.teffe))[0]
        um=np.intersect1d(um,ut)
        print('teff',len(um))

    if (input.lum > -99.0):
        ut=np.where((model['lum'] > input.lum-sig*input.lume) & \
        (model['lum'] < input.lum+sig*input.lume))[0]
        um=np.intersect1d(um,ut)
        print('lum',len(um))

    if (input.dnu > 0.0):
        model_dnu = dnusun*model['fdnu']*np.sqrt(10**model['rho'])
        ut = np.where(
            (model_dnu > input.dnu - sig*input.dnue)
            & (model_dnu < input.dnu + sig*input.dnue)
        )
        ut = ut[0]
        um = np.intersect1d(um, ut)
        print('dnu', len(um))

    if (input.numax > 0.0):
        model_numax = (numaxsun
                       * (10**model['logg']/gsun)
                       * (model['teff']/teffsun)**(-0.5))
        ut = np.where(
            (model_numax > input.numax - sig*input.numaxe)
            & (model_numax < input.numax + sig*input.numaxe)
        )
        ut = ut[0]
        um = np.intersect1d(um, ut)
        print('numax', len(um))

    if (input.logg > -99.0):
        ut = np.where(
            (model['logg'] > input.logg - sig*input.logge)
            & (model['logg'] < input.logg + sig*input.logge)
        )
        ut = ut[0]
        um = np.intersect1d(um, ut)

    if (input.feh > -99.0):
        ut = np.where(
            (model['feh'] > input.feh - sig*input.fehe)
            & (model['feh'] < input.feh + sig*input.fehe)
        )
        ut = ut[0]
        um = np.intersect1d(um, ut)
        print('feh', len(um))

    print('number of models used within non-phot obsconstraints:', len(um))

    # bail if there are not enough good models
    if (len(um) < 10):
        return result

    # add reddening
    if (redmap > -99.0):

        # if no reddening map is provided, add Av as a new variable
        # and fit for it
        if (isinstance(dustmodel,pd.DataFrame) == False):
            avs = np.arange(-0.3,1.0,0.01)

            # user-specified reddening
            #if (useav > -99.0):
            #    avs = np.zeros(1) + useav
            mod = reddening(model, um, avs, extfactors)

        # otherwise, just redden each model according to the provided map
        else:
            mod = reddening_map(
                model, model_mabs, redmap, dustmodel, um, input, extfactors, band
            )

        # photometry to use for distance
        mod_mabs = mod[band]

        um = np.arange(0,len(mod['teff']),1)

        mod['dis'] = 10**((redmap - mod_mabs + 5.0)/5.0)
        print('number of models incl reddening:',len(um))
    else:
        mod = model

    # next, another model down-select based on reddening-dependent quantities
    # only do this if no spec constraints are available
    if (mabs > -99.0):
        ut = np.where(
            (mod_mabs > mabs - sig*mabse)
            & (mod_mabs < mabs + sig*mabse)
        )
        ut = ut[0]
        um = np.intersect1d(um, ut)

    if (input.teff == -99.0):
        if ((input.bmag > -99.0) & (input.vmag > -99.0)):
            ut=np.where(
                (mod['bmag'] - mod['vmag'] > bvcol - sig*bvcole)
                & (mod['bmag'] - mod['vmag'] < bvcol + sig*bvcole))
            ut = ut[0]
            um = np.intersect1d(um,ut)

        if ((input.btmag > -99.0) & (input.vtmag > -99.0)):
            ut=np.where(
                (mod['btmag'] - mod['vtmag'] > bvtcol - sig*bvtcole)
                & (mod['btmag'] - mod['vtmag'] < bvtcol + sig*bvtcole))
            ut = ut[0]
            um = np.intersect1d(um,ut)

        if ((input.gmag > -99.0) & (input.rmag > -99.0)):
            ut = np.where(
                (mod['gmag'] - mod['rmag'] > grcol-sig*grcole)
                & (mod['gmag'] - mod['rmag'] < grcol+sig*grcole))
            ut = ut[0]
            um = np.intersect1d(um, ut)

        if ((input.rmag > -99.0) & (input.imag > -99.0)):
            ut = np.where(
                (mod['rmag'] - mod['imag'] > ricol - sig*ricole)
                & (mod['rmag'] - mod['imag'] < ricol + sig*ricole)
            )
            ut = ut[0]
            um = np.intersect1d(um,ut)

        if ((input.imag > -99.0) & (input.zmag > -99.0)):
            ut = np.where(
                (mod['imag'] - mod['zmag'] > izcol - sig*izcole)
                & (mod['imag'] - mod['zmag'] < izcol + sig*izcole)
            )
            ut = ut[0]
            um = np.intersect1d(um, ut)

        if ((input.gmag > -99.0) & (input.imag > -99.0)):
            ut = np.where(
                (mod['gmag'] - mod['imag'] > gicol-sig*gicole)
                & (mod['gmag'] - mod['imag'] < gicol+sig*gicole))
            ut = ut[0]
            um = np.intersect1d(um, ut)

        if ((input.rmag > -99.0) & (input.zmag > -99.0)):
            ut = np.where(
                (mod['rmag'] - mod['zmag'] > rzcol-sig*rzcole)
                & (mod['rmag'] - mod['zmag'] < rzcol+sig*rzcole))
            ut = ut[0]
            um = np.intersect1d(um, ut)

        if ((input.gmag > -99.0) & (input.zmag > -99.0)):
            ut = np.where(
                (mod['gmag'] - mod['zmag'] > gzcol-sig*gzcole)
                & (mod['gmag'] - mod['zmag'] < gzcol+sig*gzcole))
            ut = ut[0]
            um = np.intersect1d(um, ut)

        if ((input.jmag > -99.0) & (input.hmag > -99.0)):
            ut = np.where(
                (mod['jmag'] - mod['hmag'] > jhcol - sig*jhcole)
                & (mod['jmag'] - mod['hmag'] < jhcol + sig*jhcole)
            )
            ut = ut[0]
            um = np.intersect1d(um, ut)

        if ((input.hmag > -99.0) & (input.kmag > -99.0)):
            ut = np.where(
                (mod['hmag'] - mod['kmag'] > hkcol - sig*hkcole)
                & (mod['hmag'] - mod['kmag'] < hkcol + sig*hkcole))
            ut = ut[0]
            um = np.intersect1d(um,ut)

        if ((input.jmag > -99.0) & (input.kmag > -99.0)):
            ut = np.where(
                (mod['jmag'] - mod['kmag'] > jkcol - sig*jkcole)
                & (mod['jmag'] - mod['kmag'] < jkcol + sig*jkcole)
            )
            ut = ut[0]
            um = np.intersect1d(um, ut)

        if ((input.bpmag > -99.0) & (input.gamag > -99.0)):
            ut = np.where(
                (mod['bpmag'] - mod['gamag'] > bpgacol - sig*bpgacole)
                & (mod['bpmag'] - mod['gamag'] < bpgacol + sig*bpgacole)
            )
            ut = ut[0]
            um = np.intersect1d(um, ut)

        if ((input.gamag > -99.0) & (input.rpmag > -99.0)):
            ut = np.where(
                (mod['gamag'] - mod['rpmag'] > garpcol - sig*garpcole)
                & (mod['gamag'] - mod['rpmag'] < garpcol + sig*garpcole)
            )
            ut = ut[0]
            um = np.intersect1d(um, ut)

        if ((input.bpmag > -99.0) & (input.rpmag > -99.0)):
            ut = np.where(
                (mod['bpmag'] - mod['rpmag'] > bprpcol - sig*bprpcole)
                & (mod['bpmag'] - mod['rpmag'] < bprpcol + sig*bprpcole)
            )
            ut = ut[0]
            um = np.intersect1d(um, ut)

        if ((input.vmag > -99.0) & (input.jmag > -99.0)):
            ut = np.where(
                (mod['vmag'] - mod['jmag'] > vjcol - sig*vjcole)
                & (mod['vmag'] - mod['jmag'] < vjcol + sig*vjcole)
            )
            ut = ut[0]
            um = np.intersect1d(um, ut)

        if ((input.vtmag > -99.0) & (input.jmag > -99.0)):
            ut = np.where(
                (mod['vtmag'] - mod['jmag'] > vtjcol - sig*vtjcole)
                & (mod['vtmag'] - mod['jmag'] < vtjcol + sig*vtjcole)
            )
            ut = ut[0]
            um = np.intersect1d(um, ut)

        if ((input.gmag > -99.0) & (input.jmag > -99.0)):
            ut = np.where(
                (mod['gmag'] - mod['jmag'] > gjcol - sig*gjcole)
                & (mod['gmag'] - mod['jmag'] < gjcol + sig*gjcole)
            )
            ut = ut[0]
            um = np.intersect1d(um, ut)

        if ((input.rmag > -99.0) & (input.jmag > -99.0)):
            ut = np.where(
                (mod['rmag'] - mod['jmag'] > rjcol - sig*rjcole)
                & (mod['rmag'] - mod['jmag'] < rjcol + sig*rjcole)
            )
            ut = ut[0]
            um = np.intersect1d(um, ut)

        if ((input.vmag > -99.0) & (input.kmag > -99.0)):
            ut = np.where(
                (mod['vmag'] - mod['kmag'] > vkcol - sig*vkcole)
                & (mod['vmag'] - mod['kmag'] < vkcol + sig*vkcole)
            )
            ut = ut[0]
            um = np.intersect1d(um, ut)

        if ((input.vtmag > -99.0) & (input.kmag > -99.0)):
            ut = np.where(
                (mod['vtmag'] - mod['kmag'] > vtkcol - sig*vtkcole)
                & (mod['vtmag'] - mod['kmag'] < vtkcol + sig*vtkcole)
            )
            ut = ut[0]
            um = np.intersect1d(um, ut)

        if ((input.gmag > -99.0) & (input.kmag > -99.0)):
            ut = np.where(
                (mod['gmag'] - mod['kmag'] > gkcol - sig*gkcole)
                & (mod['gmag'] - mod['kmag'] < gkcol + sig*gkcole)
            )
            ut = ut[0]
            um = np.intersect1d(um, ut)

        if ((input.rmag > -99.0) & (input.kmag > -99.0)):
            ut = np.where(
                (mod['rmag'] - mod['kmag'] > rkcol - sig*rkcole)
                & (mod['rmag'] - mod['kmag'] < rkcol + sig*rkcole)
            )
            ut = ut[0]
            um = np.intersect1d(um, ut)

        if ((input.gamag > -99.0) & (input.jmag > -99.0)):
            ut = np.where(
                (mod['gamag'] - mod['jmag'] > gajcol - sig*gajcole)
                & (mod['gamag'] - mod['jmag'] < gajcol + sig*gajcole)
            )
            ut = ut[0]
            um = np.intersect1d(um, ut)

        if ((input.gamag > -99.0) & (input.kmag > -99.0)):
            ut = np.where(
                (mod['gamag'] - mod['kmag'] > gakcol - sig*gakcole)
                & (mod['gamag'] - mod['kmag'] < gakcol + sig*gakcole)
            )
            ut = ut[0]
            um = np.intersect1d(um, ut)

    print('number of models after phot constraints:',len(um))
    print('----')

    # bail if there are not enough good models
    if (len(um) < 10):
        return result

    def gaussian(x, mu, sig):
        return np.exp(-(x-mu)**2./(2.*sig**2.))


    # likelihoods
    if ((input.gmag > -99.0) & (input.rmag > -99.0)):
        lh_gr = gaussian(grcol, mod['gmag'][um]-mod['rmag'][um], grcole)

    else:
        lh_gr = np.ones(len(um))

    if ((input.rmag > -99.0) & (input.imag > -99.0)):
        lh_ri = gaussian(ricol, mod['rmag'][um]-mod['imag'][um], ricole)

    else:
        lh_ri = np.ones(len(um))

    if ((input.imag > -99.0) & (input.zmag > -99.0)):
        lh_iz = gaussian(izcol, mod['imag'][um]-mod['zmag'][um], izcole)
    else:
        lh_iz = np.ones(len(um))

    if ((input.gmag > -99.0) & (input.imag > -99.0)):
        lh_gi = gaussian(gicol, mod['gmag'][um]-mod['imag'][um], gicole)

    else:
        lh_gi = np.ones(len(um))

    if ((input.rmag > -99.0) & (input.zmag > -99.0)):
        lh_rz = gaussian(rzcol, mod['rmag'][um]-mod['zmag'][um], rzcole)

    else:
        lh_rz = np.ones(len(um))

    if ((input.gmag > -99.0) & (input.zmag > -99.0)):
        lh_gz = gaussian(gzcol, mod['gmag'][um]-mod['zmag'][um], gzcole)

    else:
        lh_gz = np.ones(len(um))

    if ((input.jmag > -99.0) & (input.hmag > -99.0)):
        lh_jh = gaussian(jhcol, mod['jmag'][um]-mod['hmag'][um], jhcole)
    else:
        lh_jh = np.ones(len(um))

    if ((input.hmag > -99.0) & (input.kmag > -99.0)):
        lh_hk = gaussian(hkcol, mod['hmag'][um]-mod['kmag'][um], hkcole)
    else:
        lh_hk = np.ones(len(um))

    if ((input.jmag > -99.0) & (input.kmag > -99.0)):
        lh_jk = gaussian(jkcol, mod['jmag'][um]-mod['kmag'][um], jkcole)
    else:
        lh_jk = np.ones(len(um))

    if ((input.bpmag > -99.0) & (input.gamag > -99.0)):
        lh_bpga = gaussian(bpgacol, mod['bpmag'][um]-mod['gamag'][um], bpgacole)
    else:
        lh_bpga = np.ones(len(um))

    if ((input.gamag > -99.0) & (input.rpmag > -99.0)):
        lh_garp = gaussian(garpcol, mod['gamag'][um]-mod['rpmag'][um], garpcole)
    else:
        lh_garp = np.ones(len(um))

    if ((input.bpmag > -99.0) & (input.rpmag > -99.0)):
        lh_bprp = gaussian(bprpcol, mod['bpmag'][um]-mod['rpmag'][um], bprpcole)
    else:
        lh_bprp = np.ones(len(um))

    if ((input.bmag > -99.0) & (input.vmag > -99.0)):
        lh_bv = gaussian(bvcol, mod['bmag'][um]-mod['vmag'][um], bvcole)

    else:
        lh_bv = np.ones(len(um))

    if ((input.btmag > -99.0) & (input.vtmag > -99.0)):
        lh_bvt = gaussian(bvtcol, mod['btmag'][um]-mod['vtmag'][um], bvtcole)

    else:
        lh_bvt = np.ones(len(um))

    if ((input.vmag > -99.0) & (input.jmag > -99.0)):
        lh_vj = gaussian(vjcol, mod['vmag'][um]-mod['jmag'][um], vjcole)
    else:
        lh_vj = np.ones(len(um))

    if ((input.vtmag > -99.0) & (input.jmag > -99.0)):
        lh_vtj = gaussian(vtjcol, mod['vtmag'][um]-mod['jmag'][um], vtjcole)
    else:
        lh_vtj = np.ones(len(um))

    if ((input.gmag > -99.0) & (input.jmag > -99.0)):
        lh_gj = gaussian(gjcol, mod['gmag'][um]-mod['jmag'][um], gjcole)
    else:
        lh_gj = np.ones(len(um))

    if ((input.rmag > -99.0) & (input.jmag > -99.0)):
        lh_rj = gaussian(rjcol, mod['rmag'][um]-mod['jmag'][um], rjcole)
    else:
        lh_rj = np.ones(len(um))

    if ((input.vmag > -99.0) & (input.kmag > -99.0)):
        lh_vk = gaussian(vkcol, mod['vmag'][um]-mod['kmag'][um], vkcole)
    else:
        lh_vk = np.ones(len(um))

    if ((input.vtmag > -99.0) & (input.kmag > -99.0)):
        lh_vtk = gaussian(vtkcol, mod['vtmag'][um]-mod['kmag'][um], vtkcole)
    else:
        lh_vtk = np.ones(len(um))

    if ((input.gmag > -99.0) & (input.kmag > -99.0)):
        lh_gk = gaussian(gkcol, mod['gmag'][um]-mod['kmag'][um], gkcole)
    else:
        lh_gk = np.ones(len(um))

    if ((input.rmag > -99.0) & (input.kmag > -99.0)):
        lh_rk = gaussian(rkcol, mod['rmag'][um]-mod['kmag'][um], rkcole)
    else:
        lh_rk = np.ones(len(um))

    if ((input.gamag > -99.0) & (input.jmag > -99.0)):
        lh_gaj = gaussian(gajcol, mod['gamag'][um]-mod['jmag'][um], gajcole)
    else:
        lh_gaj = np.ones(len(um))

    if ((input.gamag > -99.0) & (input.kmag > -99.0)):
        lh_gak = gaussian(gakcol, mod['gamag'][um]-mod['kmag'][um], gakcole)
    else:
        lh_gak = np.ones(len(um))

    if (input.teff > -99):
        lh_teff = gaussian(input.teff, mod['teff'][um], input.teffe)
    else:
        lh_teff = np.ones(len(um))

    if (input.lum > -99):
        lh_lum = gaussian(input.lum, mod['lum'][um], input.lume)
    else:
        lh_lum = np.ones(len(um))

    if (input.logg > -99.0):
        lh_logg = gaussian(input.logg, mod['logg'][um], input.logge)
    else:
        lh_logg = np.ones(len(um))

    if (input.feh > -99.0):
        lh_feh = gaussian(input.feh, mod['feh'][um], input.fehe)

    else:
        lh_feh = np.ones(len(um))

    if (input.plx > -99.0):
        # Compute most likely value of absolute magnitude:
        mabsIndex = np.argmax(np.exp( (-1./(2.*input.plxe**2))*(input.plx-1./mod['dis'][um])**2))

        # Only use downselected models based on input parameters:
        downSelMagArr = mod_mabs[um]

        # Compute the likelihood of the maximum magnitude given computed errors:
        lh_mabs = gaussian(downSelMagArr[mabsIndex],mod_mabs[um],mabse)
    else:
        lh_mabs = np.ones(len(um))

    if (input.dnu > 0.):
        mod_dnu = dnusun*mod['fdnu']*np.sqrt(10**mod['rho'])
        lh_dnu = np.exp( -(input.dnu-mod_dnu[um])**2.0 / (2.0*input.dnue**2.0))
    else:
        lh_dnu = np.ones(len(um))

    if (input.numax > 0.):
        mod_numax = (numaxsun
                     * (10**mod['logg']/gsun)
                     * (mod['teff']/teffsun)**(-0.5))

        lh_numax = gaussian(input.numax,mod_numax[um],input.numaxe)
    else:
        lh_numax = np.ones(len(um))

    tlh = (lh_gr*lh_ri*lh_iz*lh_gi*lh_rz*lh_gz*lh_jh*lh_hk*lh_jk*lh_bv*lh_bvt*lh_bpga*lh_garp*lh_bprp*
        lh_vj*lh_vtj*lh_gj*lh_rj*lh_vk*lh_vtk*lh_gk*lh_rk*lh_gaj*lh_gak*
        lh_teff*lh_logg*lh_feh*lh_mabs*lh_dnu*lh_numax*lh_lum)


    # metallicity prior (only if no FeH input is given)
    if (input.feh > -99.0):
        fprior = np.ones(len(um))
    else:
        fprior = fehprior(mod['feh'][um])

    # distance prior
    if (input.plx > -99.0):
        lscale = 1350.
        dprior = ((mod['dis'][um]**2/(2.0*lscale**3.))
                  *np.exp(-mod['dis'][um]/lscale))
    else:
        dprior = np.ones(len(um))

    # isochrone prior (weights)
    tprior = mod['dage'][um]*mod['dmass'][um]*mod['dfeh'][um]

    # posterior
    prob = fprior*dprior*tprior*tlh
    prob = prob/np.sum(prob)
    if (isinstance(dustmodel,pd.DataFrame) == False):
        names = ['teff', 'logg', 'feh', 'rad', 'mass', 'rho', 'lum', 'age']
        steps = [0.001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
        fixes = [0, 1, 1, 0, 0, 1, 1, 0, 1]

        if (redmap > -99.0):
            names = [
                'teff', 'logg', 'feh', 'rad', 'mass', 'rho', 'lum', 'age',
                'avs'
            ]
            steps = [0.001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
            fixes=[0, 1, 1, 0, 0, 1, 1, 0, 1]

        if ((input.plx == -99.0) & (redmap > -99)):
            names=[
                'teff', 'logg', 'feh', 'rad', 'mass', 'rho', 'lum', 'age',
                'avs', 'dis'
            ]
            steps=[0.001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
            fixes=[0, 1, 1, 0, 0, 1, 1, 0, 1, 0]

        #if ((input.plx == -99.0) & (map > -99) & (useav > -99.0)):
        #    names=['teff','logg','feh','rad','mass','rho','lum','age','dis']
        #    steps=[0.001,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]
        #    fixes=[0,1,1,0,0,1,1,0,0]

    else:
        #names=['teff','logg','feh','rad','mass','rho','lum','age']
        #steps=[0.001,0.01,0.01,0.01,0.01,0.01,0.01,0.01]
        #fixes=[0,1,1,0,0,1,1,0,1]
        #if (input.plx == -99.0):

        avstep=((np.max(mod['avs'][um])-np.min(mod['avs'][um]))/10.)
        #pdb.set_trace()

        names = [
            'teff', 'logg', 'feh', 'rad', 'mass', 'rho', 'lum', 'age', 'avs',
           'dis'
        ]
        steps=[0.001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, avstep, 0.01]
        fixes=[0, 1, 1, 0, 0, 1, 1, 0, 1, 0]


    # Provision figure
    if plot:
        plotinit()

    ix = 1
    iy = 2
    npar = len(names)

    for j in range(0,npar):
        if fnmatch.fnmatch(names[j],'*lum*'):
            lum=np.log10((mod['rad'][um]**2. * (mod['teff'][um]/5772.)**4.))
            x, y, res, err1, err2 = getpdf(
                lum, prob, name=names[j], step=steps[j], fixed=fixes[j],
                dustmodel=dustmodel)
        else:
            if (len(np.unique(mod[names[j]][um])) > 1):
                x, y, res, err1, err2 = getpdf(
                    mod[names[j]][um], prob, name=names[j], step=steps[j],
                    fixed=fixes[j],dustmodel=dustmodel
                )
            elif ((len(np.unique(mod[names[j]][um])) == 1) and (names[j] == 'avs')):
                res = mod[names[j]][um[0]]
                err1 = 0.0
                err2 = 0.0
                x = res
                y = 1.0
            else:
                res = 0.0
                err1 = 0.0
                err2 = 0.0

        print(names[j], res, err1, err2)
        setattr(result, names[j], res)
        setattr(result, names[j]+'ep', err1)
        setattr(result, names[j]+'em', err2)
        setattr(result, names[j]+'px', x)
        setattr(result, names[j]+'py', y)

        # Plot individual posteriors
        if plot:
            plotposterior(x, y, res, err1, err2, names, j, ix, iy)
            ix += 2
            iy += 2


    # calculate posteriors for a secondary with a given delta_mag, assuming it has the same
    # distance, age, and metallicity. to do this we'll interpolate the physical properties
    # of the secondary given a delta_mag, and assign it the same posterior probabilities
    # same procedure as used in Kraus+ 16
    if (input.dmag > -99.):
        print(' ')
        print('calculating properties for secondary ...')

        delta_k=input.dmag
        delta_k_err=input.dmage
        print('using dmag=',delta_k,'+/-',delta_k_err,' in ',band)

        # interpolate across constant age and metallicity
        feh_un=np.unique(mod['feh_init'][um])
        age_un=np.unique(mod['age'][um])

        #adding in the contrast error without sampling is tricky, because that uncertainty
        # is not present in the primary posterior; instead, calculate the secondary
        # posteriors 3 times for +/- contrast errors, and then add those in quadrature
        # *explicitly assumes that the contrast errors are gaussian*
        mds=[delta_k+delta_k_err,delta_k,delta_k-delta_k_err]

        # the new model quantities for the secondary
        mod_sec=np.zeros((5,3,len(prob)))

        # Now reduce model to only those that match metallicity, age, and mass (must be less than max primary mass) conditions:
        ufeh = np.in1d(model['feh_init'],feh_un) # Must match all potential primary initial metallicities
        uage = np.in1d(model['age'],age_un) # Must match all potential primary ages
        umass = np.where(model['mass'] < np.max(mod['mass'][um]))[0] # Must be less than max primary mass
        ufa = np.where((ufeh == True) & (uage == True))[0] # Find intersection of age and feh
        ufam = np.intersect1d(umass,ufa) # Find intersection of mass and ufa
        modelMin = dict((k, model[k][ufam]) for k in model.keys()) # Define minimal model grid

        # insanely inefficient triple loop follows
        for s in range(0,len(mds)):
            for r in range(0,len(feh_un)):
                for k in range (0,len(age_un)):
                    # NB the next line uses model instead of mod, since the interpolation needs
                    # the full model grid rather than the pre-selected models returned by the
                    # reddening routine (which excludes secondary solutions). This may screw
                    # things up when trying to constrain reddening (i.e. dust="none")
                    ux=np.where((modelMin['feh_init'] == feh_un[r]) & (modelMin['age'] == age_un[k]))[0]
                    ux2=np.where((mod['feh_init'][um] == feh_un[r]) & (mod['age'][um] == age_un[k]))[0]
                    sr=np.argsort(modelMin[band][ux])
                    if ((len(ux) == 0) | (len(ux2) == 0)):
                        continue
                    mod_sec[0,s,ux2]=np.interp(mod[band][um[ux2]]+mds[s],modelMin[band][ux[sr]],modelMin['teff'][ux[sr]])
                    mod_sec[1,s,ux2]=np.interp(mod[band][um[ux2]]+mds[s],modelMin[band][ux[sr]],modelMin['logg'][ux[sr]])
                    mod_sec[2,s,ux2]=np.interp(mod[band][um[ux2]]+mds[s],modelMin[band][ux[sr]],modelMin['rad'][ux[sr]])
                    mod_sec[3,s,ux2]=np.interp(mod[band][um[ux2]]+mds[s],modelMin[band][ux[sr]],modelMin['mass'][ux[sr]])
                    mod_sec[4,s,ux2]=np.interp(mod[band][um[ux2]]+mds[s],modelMin[band][ux[sr]],modelMin['rho'][ux[sr]])

        # now get PDFs across all delta mags, add errors in quadrature
        names = ['teff', 'logg', 'rad', 'mass', 'rho']
        steps=[0.001, 0.01, 0.01, 0.01, 0.01]
        fixes=[0, 1, 0, 0, 1]

        ix = 1
        iy = 2
        npar = len(names)
        for j in range(0,5):
            x, y, res_1, err1_1, err2_1 = getpdf(mod_sec[j,0,:], prob, name=names[j], step=steps[j], fixed=fixes[j],dustmodel=dustmodel)
            xo, yo, res_2, err1_2, err2_2 = getpdf(mod_sec[j,1,:], prob, name=names[j], step=steps[j], fixed=fixes[j],dustmodel=dustmodel)
            x, y, res_3, err1_3, err2_3 = getpdf(mod_sec[j,2,:], prob, name=names[j], step=steps[j], fixed=fixes[j],dustmodel=dustmodel)

            finerr1=np.sqrt(err1_2**2 + (np.abs(res_2-res_1))**2)
            finerr2=np.sqrt(err2_2**2 + (np.abs(res_2-res_3))**2)

            print(names[j], res_2, finerr1, finerr2)
            setattr(result, names[j]+'sec', res_2)
            setattr(result, names[j]+'sec'+'ep', finerr1)
            setattr(result, names[j]+'sec'+'em', finerr2)
            setattr(result, names[j]+'sec'+'px', x)
            setattr(result, names[j]+'sec'+'py', y)

            # Plot individual posteriors
            if plot:
                plotposterior_sec(xo,yo, res_2, finerr1, finerr2, names, j, ix, iy)
                ix += 2
                iy += 2


    # Plot HR diagrams
    if plot:
        plothrd(model,mod,um,input,mabs,mabse,ix,iy)

    return result

# add extinction as a model parameter
def reddening(model,um,avs,extfactors):

    model2=dict((k, model[k][um]) for k in model.keys())
    nmodels=len(model2['teff'])*len(avs)

    keys = [
        'dage', 'dmass', 'dfeh', 'teff', 'logg', 'feh', 'rad', 'mass',
        'rho', 'age', 'gmag', 'rmag', 'imag', 'zmag', 'jmag', 'hmag',
        'bmag', 'vmag', 'btmag','vtmag', 'bpmag', 'gamag', 'rpmag',
	    'dis', 'kmag', 'avs', 'fdnu', 'feh_init'
    ]

    dtype = [(key, float) for key in keys]
    model3 = np.zeros(nmodels,dtype=dtype)

    start=0
    end=len(um)

    #print start,end
    for i in range(0,len(avs)):
        ix = np.arange(start,end,1)

        # NB: in reality, the model mags should also be Av-dependent;
        # hopefully a small effect!
        for c in 'b v g r i z j h k bt vt bp ga rp'.split():
            cmag = c + 'mag'
            ac = 'a' + c
            av = extfactors['av']
            model3[cmag][ix] = model2[cmag] + avs[i]*extfactors[ac]/av

        keys = 'teff logg feh rad mass rho age feh_init dfeh dmass dage fdnu'.split()
        for key in keys:
            model3[key][ix]=model2[key]

        model3['avs'][ix] = avs[i]
        start = start + len(um)
        end = end + len(um)
        #print(i)

    return model3

# redden model given a reddening map
def reddening_map(model, model_mabs, redmap, dustmodel, um, input, extfactors,
                  band):

    if (len(band) == 4):
        bd = band[0:1]
    else:
        bd = band[0:2]

    equ = ephem.Equatorial(
        input.ra*np.pi/180.0, input.dec*np.pi/180.0, epoch=ephem.J2000
    )
    gal = ephem.Galactic(equ)
    lon_deg = gal.lon*180./np.pi
    lat_deg = gal.lat*180./np.pi

    # zero-reddening distance
    dis = 10**((redmap-model_mabs[um]+5)/5.)

    # iterate distance and map a few times
    for i in range(0,1):
        xp = np.concatenate(
            ([0.0],np.array(dustmodel.columns[2:].str[3:],dtype='float'))
        )
        fp = np.concatenate(([0.0],np.array(dustmodel.iloc[0][2:])))
        ebvs = np.interp(x=dis, xp=xp, fp = fp)
        ext_band = extfactors['a'+bd]*ebvs
        dis=10**((redmap-ext_band-model_mabs[um]+5)/5.)



    # if no models have been pre-selected (i.e. input is
    # photometry+parallax only), redden all models
    if (len(um) == len(model['teff'])):
        model3 = copy.deepcopy(model)

        for c in 'b v g r i z j h k bt vt bp ga rp'.split():
            cmag = c + 'mag'
            ac = 'a' + c
            av = extfactors['av']
            model3[cmag] = model[cmag] + extfactors[ac] * ebvs

        model3['dis'] = dis
        model3['avs'] = extfactors['av']*ebvs
        #pdb.set_trace()

    # if models have been pre-selected, extract and only redden those
    else:
        model2 = dict((k, model[k][um]) for k in model.keys())
        nmodels = len(model2['teff'])
        keys = [
            'dage', 'dmass', 'dfeh', 'teff', 'logg', 'feh', 'rad', 'mass',
            'rho', 'age', 'gmag', 'rmag', 'imag', 'zmag', 'jmag', 'hmag',
            'bmag', 'vmag', 'btmag','vtmag', 'bpmag', 'gamag', 'rpmag',
	        'dis', 'kmag', 'avs', 'fdnu', 'feh_init'
        ]

        dtype = [(key, float) for key in keys]
        model3 = np.zeros(nmodels,dtype=dtype)
        for c in 'b v g r i z j h k bt vt bp ga rp'.split():
            cmag = c + 'mag'
            ac = 'a' + c
            av = extfactors['av']
            model3[cmag] = model2[cmag] + extfactors[ac] * ebvs

        model3['dis']=dis
        model3['avs']=extfactors['av']*ebvs
        keys = 'teff logg feh rad mass rho age feh_init dfeh dmass dage fdnu'.split()
        for key in keys:
            model3[key] = model2[key]

    return model3

########################### M-dwarf error computation and gK to 2% teff uncertainty computation:
def compute_extra_MK_error(abskmag):
    massPoly = np.array([-1.218087354981032275e-04,3.202749540513295540e-03,
-2.649332720970200630e-02,5.491458806424324990e-02,6.102330369026183476e-02,
6.122397810371335014e-01])

    massPolyDeriv = np.array([-6.090436774905161376e-04,1.281099816205318216e-02,
-7.947998162910602238e-02,1.098291761284864998e-01,6.102330369026183476e-02])

    kmagExtraErr = abs(0.021*np.polyval(massPoly,abskmag)/np.polyval(massPolyDeriv,abskmag))

    return kmagExtraErr

def compute_extra_gk_color_error(gk):
    teffPoly = np.array([5.838899127633915245e-06,-4.579640759410575821e-04,
1.591988911769273360e-02,-3.229622768514631148e-01,4.234782988549875782e+00,
-3.752421323678526477e+01,2.279521336429464498e+02,-9.419602441779162518e+02,
2.570487048729761227e+03,-4.396474893847861495e+03,4.553858427460818348e+03,
-4.123317864249115701e+03,9.028586421378711748e+03])

    teffPolyDeriv = np.array([7.006678953160697955e-05,-5.037604835351633566e-03,
1.591988911769273429e-01,-2.906660491663167978e+00,3.387826390839900625e+01,
-2.626694926574968463e+02,1.367712801857678642e+03,-4.709801220889581600e+03,
1.028194819491904491e+04,-1.318942468154358357e+04,9.107716854921636696e+03,
-4.123317864249115701e+03])

    gkExtraColorErr = abs(0.02*np.polyval(teffPoly,gk)/np.polyval(teffPolyDeriv,gk))

    return gkExtraColorErr

######################################### misc stuff

# calculate parallax for each model
def redden(redmap, mabs, gl, gb, dust):
    logd = (redmap-mabs+5.)/5.
    newd = logd

    for i in range(0,1):
        cur = 10**newd
        ebv = dust(gl,gb,cur/1000.)
        av = ebv*3.1
        aj = av*1.2348743
        newd = (redmap-mabs-aj+5.)/5.

    s_newd = np.sqrt( (0.2*0.01)**2 + (0.2*0.03)**2 + (0.2*0.02)**2 )
    plx=1./(10**newd)
    s_plx=10**(-newd)*np.log(10)*s_newd
    pdb.set_trace()

    return 1./(10**newd)

def readinput(input):
    input = ascii.read('input.txt')
    ra = input['col1'][0]
    dec = input['col2'][0]
    bmag = input['col1'][1]
    bmage = input['col2'][1]
    vmag = input['col1'][2]
    vmage = input['col2'][2]
    gmag = input['col1'][3]
    gmage = input['col2'][3]
    rmag = input['col1'][4]
    rmage = input['col2'][4]
    imag = input['col1'][5]
    image = input['col2'][5]
    zmag = input['col1'][6]
    zmage = input['col2'][6]
    jmag = input['col1'][7]
    jmage = input['col2'][7]
    hmag = input['col1'][8]
    hmage = input['col2'][8]
    kmag = input['col1'][9]
    kmage = input['col2'][9]
    plx = input['col1'][10]
    plxe = input['col2'][10]
    teff = input['col1'][11]
    teffe = input['col2'][11]
    logg = input['col1'][12]
    logge = input['col2'][12]
    feh = input['col1'][13]
    fehe = input['col2'][13]
    out = (
        ra, dec, bmag, bmage, vmag, vmage, gmag, gmage, rmag, rmage,
        imag, image, zmag, zmage, jmag, jmage, hmag, hmage, kmag, kmage,
        plx, plxe, teff, teffe, logg, logge, feh, fehe
    )

    return out
