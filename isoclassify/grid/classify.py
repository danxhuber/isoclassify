import copy
import time

import ephem
import pandas as pd
import numpy as np
from astropy.io import ascii

from pdf import *  # part of isoclassify package (to do make explicit import) 
from priors import * # part of isoclassify package (to do make explicit import) 
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
        
        self.bmag = -99.0
        self.bmage = -99.0
        self.vmag = -99.0
        self.vmage = -99.0

        self.btmag = -99.0
        self.btmage = -99.0
        self.vtmag = -99.0
        self.vtmage = -99.0

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
        
    def addplx(self,value,sigma):
        self.plx = value
        self.plxe = sigma
        
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


def classify(input, model, dustmodel=0, plot=1, useav=-99.0, ext=-99.0):
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
    teffsun = 5777.0

    # bolometric correction error; kinda needs to be motivated better ...
    bcerr = 0.03

    ## extinction coefficients
    extfactors = ext
    
    ## class containing output results
    result = resdata()

    # calculate colors + errors
    bvcol = input.bmag - input.vmag
    bvtcol = input.btmag - input.vtmag
    grcol = input.gmag - input.rmag
    ricol = input.rmag - input.imag 
    izcol = input.imag - input.zmag
    jhcol = input.jmag - input.hmag
    hkcol = input.hmag - input.kmag
    bvcole = np.sqrt(input.bmage**2 + input.vmage**2)
    bvtcole = np.sqrt(input.btmage**2 + input.vtmage**2)
    grcole = np.sqrt(input.gmage**2 + input.rmage**2)
    ricole = np.sqrt(input.rmage**2 + input.image**2)
    izcole = np.sqrt(input.image**2 + input.zmage**2)
    jhcole = np.sqrt(input.jmage**2 + input.hmage**2)
    hkcole = np.sqrt(input.hmage**2 + input.kmage**2)

    # determine apparent mag to use for distance estimation. K>J>g>Vt>V
    map = -99.0
    if (input.vmag > -99.0):
        map = input.vmag
        mape = input.vmage
        band = 'v'
        model_mabs = model['vmag']

    if (input.vtmag > -99.0):
        map = input.vtmag
        mape = input.vtmage
        model_mabs = model['vtmag']
        band = 'vt'

    if (input.gmag > -99.0):
        map = input.gmag
        mape = input.gmage
        model_mabs = model['gmag']   
        band = 'g'
	
    if (input.jmag > -99.0):
        map = input.jmag
        mape = input.jmage
        model_mabs = model['jmag']
        band = 'j'

    if (input.kmag > -99.0):
        map = input.kmag
        mape = input.kmage
        model_mabs = model['kmag']
        band = 'k'
        
    # absolute magnitude
    if (input.plx > -99.0):
        mabs = -5.0 * np.log10(1.0 / input.plx) + map + 5.0
        mabse = np.sqrt(
            (-5.0 / (input.plx * np.log(10)))**2 * input.plxe**2 
            + mape**2 + bcerr**2
        )
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
        print 'teff',len(um)

    if (input.dnu > 0.0):
        model_dnu = dnusun*model['fdnu']*np.sqrt(10**model['rho'])
        ut = np.where(
            (model_dnu > input.dnu - sig*input.dnue)  
            & (model_dnu < input.dnu + sig*input.dnue)
        )
        ut = ut[0]
        um = np.intersect1d(um, ut)
        print 'dnu', len(um)

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
        print 'numax', len(um)
        
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
        print 'feh', len(um)
               
    print 'number of models used within non-phot obsconstraints:', len(um)

    # bail if there are not enough good models
    if (len(um) < 10):
        return result

    # add reddening
    if (map > -99.0):

        # if no reddening map is provided, add Av as a new variable
        # and fit for it
        if (isinstance(dustmodel,pd.DataFrame) == False):
            avs = np.arange(0.0,5.0,0.1)
            
            # user-specified reddening
            #if (useav > -99.0):
            #    avs = np.zeros(1) + useav
                
            mod = reddening(model, um, avs, extfactors)

        # otherwise, just redden each model according to the provided map
        else:
            mod = reddening_map(
                model, model_mabs, map, dustmodel, um, input, extfactors, band
            )

        # photometry to use for distance
        if (input.vmag > -99.0):
            mod_mabs = mod['vmag']

        if (input.vtmag > -99.0):
            mod_mabs = mod['vtmag']

        if (input.gmag > -99.0):
            mod_mabs = mod['gmag']

        if (input.jmag > -99.0):
            mod_mabs = mod['jmag']

        if (input.kmag > -99.0):
            mod_mabs = mod['kmag']

        um = np.arange(0,len(mod['teff']),1)

        mod['dis'] = 10**((map - mod_mabs + 5.0)/5.0)
        print 'number of models incl reddening:',len(um)
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

    print 'number of models after phot constraints:',len(um)
    print '----'

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
   
    if ((input.jmag > -99.0) & (input.hmag > -99.0)):
        lh_jh = gaussian(jhcol, mod['jmag'][um]-mod['hmag'][um], jhcole)
    else:
        lh_jh = np.ones(len(um))

    if ((input.hmag > -99.0) & (input.kmag > -99.0)):
        lh_hk = gaussian(hkcol, mod['hmag'][um]-mod['kmag'][um], hkcole)
    else:
        lh_hk = np.ones(len(um))   

    if ((input.bmag > -99.0) & (input.vmag > -99.0)):
        lh_bv = gaussian(bvcol, mod['bmag'][um]-mod['vmag'][um], bvcole)

    else:
        lh_bv = np.ones(len(um))  

    if ((input.btmag > -99.0) & (input.vtmag > -99.0)):
        lh_bvt = gaussian(bvtcol, mod['btmag'][um]-mod['vtmag'][um], bvtcole)

    else:
        lh_bvt = np.ones(len(um))   

    if (input.teff > -99):
        lh_teff = gaussian(input.teff, mod['teff'][um], input.teffe)
    else:
        lh_teff = np.ones(len(um))

    if (input.logg > -99.0):
        lh_logg = gaussian(input.logg, mod['logg'][um], input.logge)
    else:
        lh_logg = np.ones(len(um))

    if (input.feh > -99.0):
        lh_feh = gaussian(input.feh, mod['feh'][um], input.fehe)

    else:
        lh_feh = np.ones(len(um))

    if (input.plx > -99.0):
        lh_mabs = np.exp( (-1./(2.*input.plxe**2))*(input.plx-1./mod['dis'][um])**2)

        #if (input.plxe/input.plx < 0.1):
        #    lh_mabs = np.exp( -(mabs-mod_mabs[um])**2. / (2.*mabse**2.))
        #else:
        #    dv=mod_mabs[um]-mabs
        #    lh_mabs = 10**(0.2*dv) * np.exp(-((10**(0.2*dv)-1.)**2)/(2.*(input.plxe/input.plx)**2))
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

    tlh = (lh_gr*lh_ri*lh_iz*lh_jh*lh_hk*lh_bv*lh_bvt*lh_teff*lh_logg*lh_feh
           *lh_mabs*lh_dnu*lh_numax)
        
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
        
        if (map > -99.0):
            names = [
                'teff', 'logg', 'feh', 'rad', 'mass', 'rho', 'lum', 'age',
                'avs'
            ]
            steps = [0.001, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
            fixes=[0, 1, 1, 0, 0, 1, 1, 0, 1]

        if ((input.plx == -99.0) & (map > -99)):
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
            lum=np.log10((mod['rad'][um]**2. * (mod['teff'][um]/5777.)**4.))
            x, y, res, err1, err2 = getpdf(
                lum, prob, name=names[j], step=steps[j], fixed=fixes[j],
                dustmodel=dustmodel)
        else:
            if (len(np.unique(mod[names[j]][um])) > 1):
                x, y, res, err1, err2 = getpdf(
                    mod[names[j]][um], prob, name=names[j], step=steps[j],
                    fixed=fixes[j],dustmodel=dustmodel
                )
            else:
                res = 0.0
                err1 = 0.0
                err2 = 0.0

        print names[j], res, err1, err2
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
    
    # Plot HR diagrams
    if plot:
        plothrd(model,input,mabs,mabse,ix,iy)
        
    # angular diameter
    rsun=6.9599e10
    pc=3.086e18
    rad2mas=206265000.
    theta = ((mod['rad']*rsun)/(mod['dis']*pc))*rad2mas*2.        
        
    x, y, res, err1, err2 = getpdf(
            theta[um], prob, name='angdia', step=0.001,
            fixed=1,dustmodel=dustmodel)
                
    print 'angdia', res, err1, err2
    plt.figure(3)
    plt.subplot(2,1,1)
    plt.plot(x,np.cumsum(y))
    plt.plot([res,res],[0,1],'r')
    plt.plot([res+err1,res+err1],[0,1],'--r')
    plt.plot([res-err2,res-err2],[0,1],'--r')
    plt.subplot(2,1,2)
    plt.plot(x,y)
    plt.plot([res,res],[0,1],'r')
    plt.plot([res+err1,res+err1],[0,1],'--r')
    plt.plot([res-err2,res-err2],[0,1],'--r')
    plt.ylim([0,np.max(y)+np.max(y)*0.1])
    
    #pdb.set_trace()

    return result
            
# add extinction as a model parameter
def reddening(model,um,avs,extfactors):

    model2=dict((k, model[k][um]) for k in model.keys())
    nmodels=len(model2['teff'])*len(avs)
    #pdb.set_trace()

    keys = [
        'dage', 'dmass', 'dfeh', 'teff', 'logg', 'feh', 'rad', 'mass', 
        'rho', 'age', 'gmag', 'rmag', 'imag', 'zmag', 'jmag', 'hmag', 
        'bmag', 'vmag', 'btmag','vtmag', 'dis', 'kmag', 'avs', 'fdnu'
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
        for c in 'b v g r i z j h k bt vt'.split():
            cmag = c + 'mag'
            ac = 'a' + c
            av = extfactors['av']
            model3[cmag][ix] = model2[cmag] + avs[i]*extfactors[ac]/av

        keys = 'teff logg feh rad mass rho age dfeh dmass dage fdnu'.split()
        for key in keys:
            model3[key][ix]=model2[key]

        model3['avs'][ix] = avs[i]
        start = start + len(um)
        end = end + len(um)

    return model3

# redden model given a reddening map
def reddening_map(model, model_mabs, map, dustmodel, um, input, extfactors, 
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
    dis = 10**((map-model_mabs[um]+5)/5.)

    # iterate distance and map a few times
    for i in range(0,1):
        xp = np.concatenate(
            ([0.0],np.array(dustmodel.columns[2:].str[3:],dtype='float'))
        )
        fp = np.concatenate(([0.0],np.array(dustmodel.iloc[0][2:])))
        ebvs = np.interp(x=dis, xp=xp, fp = fp)
        ext_band = extfactors['a'+bd]*ebvs	
        dis=10**((map-ext_band-model_mabs[um]+5)/5.)

    # if no models have been pre-selected (i.e. input is
    # photometry+parallax only), redden all models
    if (len(um) == len(model['teff'])):
        model3 = copy.deepcopy(model)

        for c in 'b v g r i z j h k bt vt'.split():
            cmag = c + 'mag'
            ac = 'a' + c
            av = extfactors['av']
            model3[cmag] = model[cmag] + extfactors[ac] * ebvs

        model3['dis'] = dis
        model3['avs'] = extfactors['av']*ebvs	
	 
    # if models have been pre-selected, extract and only redden those
    else:
        model2 = dict((k, model[k][um]) for k in model.keys())
        nmodels = len(model2['teff'])
        keys = [
            'dage', 'dmass', 'dfeh', 'teff', 'logg', 'feh', 'rad', 'mass', 
            'rho', 'age', 'gmag', 'rmag', 'imag', 'zmag', 'jmag', 'hmag', 
            'bmag', 'vmag', 'btmag','vtmag', 'dis', 'kmag', 'avs', 'fdnu'
        ]
        
        dtype = [(key, float) for key in keys]
        model3 = np.zeros(nmodels,dtype=dtype)
        for c in 'b v g r i z j h k bt vt'.split():
            cmag = c + 'mag'
            ac = 'a' + c
            av = extfactors['av']
            model3[cmag] = model2[cmag] + extfactors[ac] * ebvs

        model3['dis']=dis
        model3['avs']=extfactors['av']*ebvs	
        keys = 'teff logg feh rad mass rho age dfeh dmass dage fdnu'.split()
        for key in keys:
            model3[key] = model2[key]

    return model3


######################################### misc stuff

# calculate parallax for each model
def redden(map, mabs, gl, gb, dust):
    logd = (map-mabs+5.)/5.
    newd = logd

    for i in range(0,1):
        cur = 10**newd
        ebv = dust(gl,gb,cur/1000.)
        av = ebv*3.1
        aj = av*1.2348743
        newd = (map-mabs-aj+5.)/5.

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

