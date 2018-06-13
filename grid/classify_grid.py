from pdf import *
from priors import *
from plot import *
from astropy.io import ascii
import ephem
import pdb
import numpy as np
#import mwdust
import time
import copy
import pandas as pd

class obsdata():
    def __init__(self):
        self.plx = -99.
        self.plxe = -99.
    
        self.teff = -99.
        self.teffe = -99.
        self.logg = -99.
        self.logge = -99.
        self.feh = -99.
        self.fehe = -99.
        
        self.bmag = -99.
        self.bmage = -99.
        self.vmag = -99.
        self.vmage = -99.

        self.btmag = -99.
        self.btmage = -99.
        self.vtmag = -99.
        self.vtmage = -99.

        self.gmag = -99.
        self.gmage = -99.
        self.rmag = -99.
        self.rmage = -99.
        self.imag = -99.
        self.image = -99.
        self.zmag = -99.
        self.zmage = -99.
        self.jmag = -99.
        self.jmage = -99.
        self.hmag = -99.
        self.hmage = -99.
        self.kmag = -99.
        self.kmage = -99.
        
        self.numax = -99.
        self.numaxe = -99.
        self.dnu = -99.
        self.dnue = -99.
                   
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
        self.teff = 0.
        self.teffep = 0.
        self.teffem = 0.
        self.teffpx = 0.
        self.teffpy = 0.
        self.logg = 0.
        self.loggep = 0.
        self.loggem = 0.
        self.loggpx = 0.
        self.loggpy = 0.
        self.feh = 0.
        self.fehep = 0.
        self.fehem = 0.
        self.fehpx = 0.
        self.fehpy = 0.
        self.rad = 0.
        self.radep = 0.
        self.radem = 0.
        self.radpx = 0.
        self.radpy = 0.
        self.mass = 0.
        self.massep = 0.
        self.massem = 0.
        self.masspx = 0.
        self.masspy = 0.
        self.rho = 0.
        self.rhoep = 0.
        self.rhoem = 0.
        self.rhopx = 0.
        self.rhopy = 0.
        self.lum = 0.
        self.lumep = 0.
        self.lumem = 0.
        self.lumpx = 0.
        self.lumpy = 0.
        self.age = 0.
        self.ageep = 0.
        self.ageem = 0.
        self.agepx = 0.
        self.agepy = 0.
        self.avs = 0.
        self.avsep = 0.
        self.avsem = 0.
        self.avspx = 0.
        self.avspy = 0.
        self.dis = 0.
        self.disep = 0.
        self.disem = 0.
        self.dispx = 0.
        self.dispy = 0.

class extinction():
    def __init__(self):
        self.ab=1.3454449
        self.av=1.00

        self.abt=1.3986523
        self.avt=1.0602271

        self.ag=1.2348743
        self.ar=0.88343449
        self.ai=0.68095687
        self.az=0.48308430

        self.aj=0.28814896
        self.ah=0.18152716
        self.ak=0.11505195

        self.aga=1.2348743


def classify(input,model,dustmodel=0,doplot=1,useav=-99.):

    ## constants
    gsun=27420.010
    numaxsun=3090.
    dnusun=135.1
    teffsun=5777.

    # bolometric correction error; kinda needs to be motivated better ...
    bcerr=0.03

    ## extinction coefficients
    extfactors=extinction()
    
    ## class containing output results
    result = resdata()

    # calculate colors + errors
    bvcol=input.bmag-input.vmag
    bvtcol=input.btmag-input.vtmag
    grcol=input.gmag-input.rmag
    ricol=input.rmag-input.imag 
    izcol=input.imag-input.zmag
    jhcol=input.jmag-input.hmag
    hkcol=input.hmag-input.kmag
    bvcole=np.sqrt(input.bmage**2 + input.vmage**2)
    bvtcole=np.sqrt(input.btmage**2 + input.vtmage**2)
    grcole=np.sqrt(input.gmage**2 + input.rmage**2)
    ricole=np.sqrt(input.rmage**2 + input.image**2)
    izcole=np.sqrt(input.image**2 + input.zmage**2)
    jhcole=np.sqrt(input.jmage**2 + input.hmage**2)
    hkcole=np.sqrt(input.hmage**2 + input.kmage**2)

    #bcmodel = h5py.File('bcgrid.h5', 'r')

    # determine apparent mag to use for distance estimation. K>J>g>Vt>V
    map=-99.
    if (input.vmag > -99.):
        map=input.vmag
        mape=input.vmage
        ext=extfactors.av
        model_mabs=model['vmag']

    if (input.vtmag > -99.):
        map=input.vtmag
        mape=input.vtmage
        model_mabs=model['vtmag']
	ext=extfactors.avt

    if (input.gmag > -99.):
        map=input.gmag
        mape=input.gmage
        model_mabs=model['gmag']   
        ext=extfactors.ag
	
    if (input.jmag > -99.):
        map=input.jmag
        mape=input.jmage
        model_mabs=model['jmag']
        ext=extfactors.aj

    if (input.kmag > -99.):
        map=input.kmag
        mape=input.kmage
        model_mabs=model['kmag']
        ext=extfactors.ak
        
    # absolute magnitude
    if (input.plx > -99.):
        mabs = -5.*np.log10(1./input.plx)+map+5.
        mabse = np.sqrt((-5./(input.plx*np.log(10)))**2*input.plxe**2+mape**2+bcerr**2.)
    else:
        mabs=-99.
        mabse=-99.

    # pre-select model grid; first only using reddening-independent quantities
    sig=4.
    um=np.arange(0,len(model['teff']),1)
        
    if (input.teff > -99.):
        ut=np.where((model['teff'] > input.teff-sig*input.teffe) & \
        (model['teff'] < input.teff+sig*input.teffe))[0]
        um=np.intersect1d(um,ut)
        print 'teff',len(um)

    if (input.dnu > 0.):
        model_dnu = dnusun*model['fdnu']*np.sqrt(10**model['rho'])
        ut=np.where((model_dnu > input.dnu-sig*input.dnue) & \
        (model_dnu < input.dnu+sig*input.dnue))[0]
        um=np.intersect1d(um,ut)
        print 'dnu',len(um)

    if (input.numax > 0.):
        model_numax = numaxsun*(10**model['logg']/gsun)*(model['teff']/teffsun)**(-0.5)
        ut=np.where((model_numax > input.numax-sig*input.numaxe) & \
        (model_numax < input.numax+sig*input.numaxe))[0]
        um=np.intersect1d(um,ut)
        print 'numax',len(um)
        
    if (input.logg > -99.):
        ut=np.where((model['logg'] > input.logg-sig*input.logge) & \
        (model['logg'] < input.logg+sig*input.logge))[0]
        um=np.intersect1d(um,ut)
        
    if (input.feh > -99.):
        ut=np.where((model['feh'] > input.feh-sig*input.fehe) & \
        (model['feh'] < input.feh+sig*input.fehe))[0]
        um=np.intersect1d(um,ut)
        print 'feh',len(um)
               
    print 'number of models used within non-phot obsconstraints:',len(um)
    # bail if there are not enough good models
    if (len(um) < 10):
        return result

    # add reddening
    if (map > -99.):

        # if no reddening map is provided, add Av as a new variable and fit for it
        if (isinstance(dustmodel,pd.DataFrame) == False):
            #avs = np.linspace(-0.1,1.0,41.)
            avs = np.arange(-0.3,1.0,0.01)
            #avs = np.arange(-0.3,1.0,0.1)
            
            # user-specified reddening
            if (useav > -99.):
                avs = np.zeros(1)+useav
                
            mod = reddening(model,um,avs,extfactors)

        # otherwise, just redden each model according to the provided map
        else:
	    mod=reddening_map(model,model_mabs,map,ext,dustmodel,um,input,extfactors)
        
        # photometry to use for distance
        if (input.vmag > -99.):
            mod_mabs=mod['vmag']
        if (input.vtmag > -99.):
            mod_mabs=mod['vtmag']
	if (input.gmag > -99.):
            mod_mabs=mod['gmag']
	if (input.jmag > -99.):
            mod_mabs=mod['jmag']
        if (input.kmag > -99.):
            mod_mabs=mod['kmag']
        um = np.arange(0,len(mod['teff']),1)
	
        mod['dis'] = 10**((map-mod_mabs+5.)/5.)
        print 'number of models incl reddening:',len(um)
    else:
        mod=model

    # next, another model down-select based on reddening-dependent quantities
    # only do this if no spec constraints are available
    
    if (mabs > -99.):
            ut=np.where((mod_mabs > mabs-sig*mabse) & (mod_mabs < mabs+sig*mabse))[0]
            um=np.intersect1d(um,ut)
    #else:
    #        um = np.arange(0,len(mod['teff']),1)
    

    if (input.teff == -99.):

        if ((input.bmag > -99.) & (input.vmag > -99.)):
            ut=np.where((model['bmag']-model['vmag'] > bvcol-sig*bvcole) & \
                        (model['bmag']-model['vmag'] < bvcol+sig*bvcole))[0]
            um=np.intersect1d(um,ut)
            pdb.set_trace()
            #print 'bv:',len(um)
        if ((input.btmag > -99.) & (input.vtmag > -99.)):
            ut=np.where((model['btmag']-model['vtmag'] > bvtcol-sig*bvtcole) & \
                    (model['btmag']-model['vtmag'] < bvtcol+sig*bvtcole))[0]
            um=np.intersect1d(um,ut)
            #print 'btvt:',len(um)
        if ((input.gmag > -99.) & (input.rmag > -99.)):
            ut=np.where((model['gmag']-model['rmag'] > grcol-sig*grcole) & \
                    (model['gmag']-model['rmag'] < grcol+sig*grcole))[0]
            um=np.intersect1d(um,ut)
            #print 'gr:',len(um)
        if ((input.rmag > -99.) & (input.imag > -99.)):
            ut=np.where((model['rmag']-model['imag'] > ricol-sig*ricole) & \
                    (model['rmag']-model['imag'] < ricol+sig*ricole))[0]
            um=np.intersect1d(um,ut)
            #print 'ri:',len(um)
        if ((input.imag > -99.) & (input.zmag > -99.)):
            ut=np.where((model['imag']-model['zmag'] > izcol-sig*izcole) & \
                    (model['imag']-model['zmag'] < izcol+sig*izcole))[0]
            um=np.intersect1d(um,ut)
            #print 'iz:',len(um)
        if ((input.jmag > -99.) & (input.hmag > -99.)):
            ut=np.where((model['jmag']-model['hmag'] > jhcol-sig*jhcole) & \
                    (model['jmag']-model['hmag'] < jhcol+sig*jhcole))[0]
            um=np.intersect1d(um,ut)
            #print 'jh:',len(um)
        if ((input.hmag > -99.) & (input.kmag > -99.)):
            ut=np.where((model['hmag']-model['kmag'] > hkcol-sig*hkcole) & \
                    (model['hmag']-model['kmag'] < hkcol+sig*hkcole))[0]
            um=np.intersect1d(um,ut)
            #print 'hk:',len(um)

        
    #pdb.set_trace()

    print 'number of models after phot constraints:',len(um)
    print '----'

    # bail if there are not enough good models
    if (len(um) < 10):
        return result
    
    # likelihoods
    if ((input.gmag > -99.) & (input.rmag > -99.)):
        lh_gr = np.exp(-(grcol-(mod['gmag'][um]-mod['rmag'][um]))**2./(2.*grcole**2.))
    else:
        lh_gr = np.ones(len(um))

    if ((input.rmag > -99.) & (input.imag > -99.)):
        lh_ri = np.exp(-(ricol-(mod['rmag'][um]-mod['imag'][um]))**2./(2.*ricole**2.))
    else:
        lh_ri = np.ones(len(um)) 
        
    if ((input.imag > -99.) & (input.zmag > -99.)):
        lh_iz = np.exp(-(izcol-(mod['imag'][um]-mod['zmag'][um]))**2./(2.*izcole**2.))
    else:
        lh_iz = np.ones(len(um)) 
   
    if ((input.jmag > -99.) & (input.hmag > -99.)):
        lh_jh = np.exp(-(jhcol-(mod['jmag'][um]-mod['hmag'][um]))**2./(2.*jhcole**2.))
    else:
        lh_jh = np.ones(len(um))

    if ((input.hmag > -99.) & (input.kmag > -99.)):
        lh_hk = np.exp(-(hkcol-(mod['hmag'][um]-mod['kmag'][um]))**2./(2.*hkcole**2.))
    else:
        lh_hk = np.ones(len(um))   

    if ((input.bmag > -99.) & (input.vmag > -99.)):
        lh_bv = np.exp(-(bvcol-(mod['bmag'][um]-mod['vmag'][um]))**2./(2.*bvcole**2.))
    else:
        lh_bv = np.ones(len(um))  

    if ((input.btmag > -99.) & (input.vtmag > -99.)):
        lh_bvt = np.exp(-(bvtcol-(mod['btmag'][um]-mod['vtmag'][um]))**2./(2.*bvtcole**2.))
    else:
        lh_bvt = np.ones(len(um))   

    if (input.teff > -99):
        lh_teff = np.exp(-(input.teff-mod['teff'][um])**2./(2.*input.teffe**2.))
    else:
        lh_teff = np.ones(len(um))

    if (input.logg > -99.):
        lh_logg = np.exp(-(input.logg-mod['logg'][um])**2./(2.*input.logge**2.))
    else:
        lh_logg = np.ones(len(um))

    if (input.feh > -99.):
        lh_feh = np.exp(-(input.feh-mod['feh'][um])**2./(2.*input.fehe**2.))
    else:
        lh_feh = np.ones(len(um))

    if (input.plx > -99.):
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
        lh_dnu = np.exp( -(input.dnu-mod_dnu[um])**2. / (2.*input.dnue**2.))
    else:
        lh_dnu = np.ones(len(um))

    if (input.numax > 0.):
        mod_numax = numaxsun*(10**mod['logg']/gsun)*(mod['teff']/teffsun)**(-0.5)
        lh_numax = np.exp( -(input.numax-mod_numax[um])**2. / (2.*input.numaxe**2.))
    else:
        lh_numax = np.ones(len(um))

    tlh = lh_gr*lh_ri*lh_iz*lh_jh*lh_hk*lh_bv*lh_bvt*lh_teff*lh_logg*lh_feh*lh_mabs*lh_dnu*lh_numax
    #pdb.set_trace()
        
    # metallicity prior (only if no FeH input is given)
    if (input.feh > -99.):
        fprior = np.ones(len(um))
    else:
        fprior = fehprior(mod['feh'][um])
    
    # distance prior
    if (input.plx > -99.):
        lscale = 1350.
        dprior = (mod['dis'][um]**2/(2.*lscale**3.))*np.exp(-mod['dis'][um]/lscale)
    else:
        dprior = np.ones(len(um))

    # isochrone prior (weights)
    tprior = mod['dage'][um]*mod['dmass'][um]*mod['dfeh'][um]

    # posterior
    prob = fprior*dprior*tprior*tlh
    prob = prob/np.sum(prob)
    #pdb.set_trace()
    
    if (isinstance(dustmodel,pd.DataFrame) == False):

        names=['teff','logg','feh','rad','mass','rho','lum','age']
        steps=[0.001,0.01,0.01,0.01,0.01,0.01,0.01,0.01]
        fixes=[0,1,1,0,0,1,1,0,1]
        #pdb.set_trace()
        
        if (map > -99.):
            names=['teff','logg','feh','rad','mass','rho','lum','age','avs']
            steps=[0.001,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]
            fixes=[0,1,1,0,0,1,1,0,1]

        if ((input.plx == -99.) & (map > -99)):
            names=['teff','logg','feh','rad','mass','rho','lum','age','avs','dis']
            steps=[0.001,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]
            fixes=[0,1,1,0,0,1,1,0,1,0]
            
        if ((input.plx == -99.) & (map > -99) & (useav > -99.)):
            names=['teff','logg','feh','rad','mass','rho','lum','age','dis']
            steps=[0.001,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]
            fixes=[0,1,1,0,0,1,1,0,0]
            
        
            
    else:
        #names=['teff','logg','feh','rad','mass','rho','lum','age']
        #steps=[0.001,0.01,0.01,0.01,0.01,0.01,0.01,0.01]
        #fixes=[0,1,1,0,0,1,1,0,1]
        #if (input.plx == -99.):

        avstep=((np.max(mod['avs'][um])-np.min(mod['avs'][um]))/10.)
        #pdb.set_trace()

        names=['teff','logg','feh','rad','mass','rho','lum','age','avs','dis']
        steps=[0.001,0.01,0.01,0.01,0.01,0.01,0.01,0.01,avstep,0.01]
        fixes=[0,1,1,0,0,1,1,0,1,0]

    #pdb.set_trace()
                           
    if doplot:
            plotinit()

    ix=1
    iy=2
    npar=len(names)
   
    for j in range(0,npar):

            if fnmatch.fnmatch(names[j],'*lum*'):
                lum=np.log10((mod['rad'][um]**2. * (mod['teff'][um]/5777.)**4.))
                x,y,res,err1,err2 = getpdf(lum,prob,name=names[j],\
                                           step=steps[j],fixed=fixes[j],dustmodel=dustmodel)
            else:
                if (len(np.unique(mod[names[j]][um])) > 1):
                    x,y,res,err1,err2 = getpdf(mod[names[j]][um],prob,name=names[j],\
                                               step=steps[j],fixed=fixes[j],dustmodel=dustmodel)
                else:
                    res=0.
                    err1=0.
                    err2=0.

            print names[j],res,err1,err2
            setattr(result, names[j], res)
            setattr(result, names[j]+'ep', err1)
            setattr(result, names[j]+'em', err2)
            setattr(result, names[j]+'px', x)
            setattr(result, names[j]+'py', y)

            if doplot:
                    #plotposterior(x,y,res,err1,err2,0.,model,model,names,j,0.,\
                    #              0.,grcol,ricol,grcole,ricole,mabs,mabse,ix,iy)
                    plotposterior(x,y,res,err1,err2,names,j,ix,iy)
                    ix=ix+2
                    iy=iy+2

    if doplot:
#            plothrd(x,y,res,err1,err2,0.,model,model,names,j,0.,0.,grcol,ricol,grcole,\
#                     ricole,0.,0.,ix,iy)
            plothrd(model,input,mabs,mabse,ix,iy)

            plotclear()

    #pdb.set_trace()
    return result
            
    #if i==300:
    #pdb.set_trace()



# add extinction as a model parameter
def reddening(model,um,avs,extfactors):

    model2=dict((k, model[k][um]) for k in model.keys())
    nmodels=len(model2['teff'])*len(avs)
    #pdb.set_trace()
    model3=np.zeros(nmodels,dtype=[('dage', float), ('dmass', float), ('dfeh', float), \
            ('teff', float), ('logg', float), ('feh', float), \
            ('rad', float), ('mass', float), ('rho', float), ('age', float), \
            ('gmag', float), ('rmag', float), ('imag', float), \
            ('zmag', float), ('jmag', float), ('hmag', float), \
            ('bmag', float), ('vmag', float), ('btmag', float), ('vtmag', float), ('dis', float), \
            ('kmag', float), ('avs', float), ('fdnu', float)])

    start=0
    end=len(um)

    #print start,end
    for i in range(0,len(avs)):
            ix=np.arange(start,end,1)

            # NB: in reality, the model mags should also be Av-dependent; hopefully a 
            # small effect!
            model3['bmag'][ix]=model2['bmag']+avs[i]*extfactors.ab
            model3['vmag'][ix]=model2['vmag']+avs[i]*extfactors.av
            model3['gmag'][ix]=model2['gmag']+avs[i]*extfactors.ag
            model3['rmag'][ix]=model2['rmag']+avs[i]*extfactors.ar
            model3['imag'][ix]=model2['imag']+avs[i]*extfactors.ai
            model3['zmag'][ix]=model2['zmag']+avs[i]*extfactors.az
            model3['jmag'][ix]=model2['jmag']+avs[i]*extfactors.aj
            model3['hmag'][ix]=model2['hmag']+avs[i]*extfactors.ah
            model3['kmag'][ix]=model2['kmag']+avs[i]*extfactors.ak
            model3['btmag'][ix]=model2['btmag']+avs[i]*extfactors.abt
            model3['vtmag'][ix]=model2['vtmag']+avs[i]*extfactors.avt

            model3['teff'][ix]=model2['teff']
            model3['logg'][ix]=model2['logg']
            model3['feh'][ix]=model2['feh']
            model3['rad'][ix]=model2['rad']
            model3['mass'][ix]=model2['mass']
            model3['rho'][ix]=model2['rho']
            model3['age'][ix]=model2['age']
            model3['dfeh'][ix]=model2['dfeh']
            model3['dmass'][ix]=model2['dmass']
            model3['dage'][ix]=model2['dage']
            model3['fdnu'][ix]=model2['fdnu']

            model3['avs'][ix]=avs[i]
            start=start+len(um)
            end=end+len(um)

    return model3


# redden model given a reddening map
def reddening_map(model,model_mabs,map,ext,dustmodel,um,input,extfactors):

    equ = ephem.Equatorial(input.ra*np.pi/180., input.dec*np.pi/180., epoch=ephem.J2000)
    gal = ephem.Galactic(equ)
    lon_deg=gal.lon*180./np.pi
    lat_deg=gal.lat*180./np.pi

    # zero-reddening distance
    dis = 10**((map-model_mabs[um]+5)/5.)

    # iterate distance and map a few times
    for i in range(0,1):
        ext_v = (3.1)*np.interp(x=dis/1000.,xp=np.concatenate(([0.0],np.array(dustmodel.columns[2:].str[3:],dtype='float'))),fp=np.concatenate(([0.0],np.array(dustmodel.iloc[0][2:]))))
        ext_band = ext_v*ext	
        dis=10**((map-ext_band-model_mabs[um]+5)/5.)
  
    
    # if no models have been pre-selected (i.e. input is photometry+parallax only), 
    # redden all models
    if (len(um) == len(model['teff'])):
	    model3 = copy.deepcopy(model)
	    model3['bmag']=model['bmag']+ext_v*extfactors.ab
	    model3['vmag']=model['vmag']+ext_v*extfactors.av
	    model3['gmag']=model['gmag']+ext_v*extfactors.ag
	    model3['rmag']=model['rmag']+ext_v*extfactors.ar
	    model3['imag']=model['imag']+ext_v*extfactors.ai
	    model3['zmag']=model['zmag']+ext_v*extfactors.az
	    #model3['d51mag']=model['d51mag']+ext_v*ab
	    model3['jmag']=model['jmag']+ext_v*extfactors.aj
	    model3['hmag']=model['hmag']+ext_v*extfactors.ah
	    model3['kmag']=model['kmag']+ext_v*extfactors.ak
	    model3['btmag']=model['btmag']+ext_v*extfactors.abt
	    model3['vtmag']=model['vtmag']+ext_v*extfactors.avt
	    model3['dis']=dis
	    model3['avs']=ext_v
	    #pdb.set_trace()
	 
    # if models have been pre-selected, extract and only redden those
    else:
	    model2=dict((k, model[k][um]) for k in model.keys())
	    nmodels=len(model2['teff'])
	    model3=np.zeros(nmodels,dtype=[('dage', float), ('dmass', float), ('dfeh', float), \
        	    ('teff', float), ('logg', float), ('feh', float), \
        	    ('rad', float), ('mass', float), ('rho', float), ('age', float), \
        	    ('gmag', float), ('rmag', float), ('imag', float), \
        	    ('zmag', float), ('jmag', float), ('hmag', float), \
        	    ('bmag', float), ('vmag', float), ('btmag', float), ('vtmag', float), ('dis', float), \
        	    ('kmag', float), ('avs', float), ('fdnu', float)])

	    model3['bmag']=model2['bmag']+ext_v*extfactors.ab
	    model3['vmag']=model2['vmag']+ext_v*extfactors.av
	    model3['gmag']=model2['gmag']+ext_v*extfactors.ag
	    model3['rmag']=model2['rmag']+ext_v*extfactors.ar
	    model3['imag']=model2['imag']+ext_v*extfactors.ai
	    model3['zmag']=model2['zmag']+ext_v*extfactors.az
	    #model3['d51mag']=model2['d51mag']+ext_v*ab
	    model3['jmag']=model2['jmag']+ext_v*extfactors.aj
	    model3['hmag']=model2['hmag']+ext_v*extfactors.ah
	    model3['kmag']=model2['kmag']+ext_v*extfactors.ak
	    model3['btmag']=model2['btmag']+ext_v*extfactors.abt
	    model3['vtmag']=model2['vtmag']+ext_v*extfactors.avt
	    model3['dis']=dis
	    model3['avs']=ext_v

	    model3['teff']=model2['teff']
	    model3['logg']=model2['logg']
	    model3['feh']=model2['feh']
	    model3['rad']=model2['rad']
	    model3['mass']=model2['mass']
	    model3['rho']=model2['rho']
	    model3['age']=model2['age']
	    model3['dfeh']=model2['dfeh']
	    model3['dmass']=model2['dmass']
	    model3['dage']=model2['dage']
	    model3['fdnu']=model2['fdnu']
            #pdb.set_trace()
       
    return model3


















######################################### misc stuff

# calculate parallax for each model
def redden(map,mabs,gl,gb,dust):

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
    input=ascii.read('input.txt')
    ra=input['col1'][0]
    dec=input['col2'][0]
    bmag=input['col1'][1]
    bmage=input['col2'][1]
    vmag=input['col1'][2]
    vmage=input['col2'][2]
    gmag=input['col1'][3]
    gmage=input['col2'][3]
    rmag=input['col1'][4]
    rmage=input['col2'][4]
    imag=input['col1'][5]
    image=input['col2'][5]
    zmag=input['col1'][6]
    zmage=input['col2'][6]
    jmag=input['col1'][7]
    jmage=input['col2'][7]
    hmag=input['col1'][8]
    hmage=input['col2'][8]
    kmag=input['col1'][9]
    kmage=input['col2'][9]
    plx=input['col1'][10]
    plxe=input['col2'][10]
    teff=input['col1'][11]
    teffe=input['col2'][11]
    logg=input['col1'][12]
    logge=input['col2'][12]
    feh=input['col1'][13]
    fehe=input['col2'][13]
    return ra,dec,bmag,bmage,vmag,vmage,gmag,gmage,rmag,rmage,imag,image,zmag,zmage,\
        jmag,jmage,hmag,hmage,kmag,kmage,plx,plxe,teff,teffe,logg,logge,feh,fehe
