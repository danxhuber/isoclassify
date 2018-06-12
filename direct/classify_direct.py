# code to calculate fundamental stellar parameters and distances using 
# a "direct method", i.e. adopting a fixed reddening map and bolometric 
# corrections

import numpy as np
import pdb
import astropy.units as units
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord

def stparas(input,dnumodel=-99,bcmodel=-99,dustmodel=-99,dnucor=-99,useav=-99,plot=-99,band='k'):

    # IAU XXIX Resolution, Mamajek et al. (2015)
    r_sun = 6.957e10   	    
    gconst = 6.67408e-8	    
    gm=1.3271244e26     	    	
    m_sun=gm/gconst
    rho_sun=m_sun/(4./3.*np.pi*r_sun**3)
    g_sun = gconst*m_sun/r_sun**2.

    # solar constants
    numaxsun = 3090.
    dnusun = 135.1
    teffsun = 5777.
    Msun = 4.74 # NB this is fixed to MESA BCs!

    # assumed uncertainty in bolometric corrections
    err_bc=0.02

    # assumed uncertainty in extinction
    err_ext=0.02

    # object containing output values
    out = resdata()

    ## extinction coefficients
    #extfactors=extinction()
    extfactors={"ab":1.3454449, "av":1.0, "abt":1.3454449, "avt":1.0602271, "ag":1.2348743, "ar":0.88343449, "ai":0.68095687, "az":0.48308430, "aj":0.28814896, "ah":0.18152716, "ak":0.11505195, "aga":1.2348743}

    if (len(band) == 4):
        bd=band[0:1]
    else:
        bd=band[0:2]

    ######################################
    # case 1: input is parallax + colors #
    ######################################

    #with h5py.File(bcmodel,'r') as h5:
    teffgrid = bcmodel['teffgrid'][:]
    logggrid = bcmodel['logggrid'][:]
    fehgrid = bcmodel['fehgrid'][:]
    avgrid = bcmodel['avgrid'][:]
    bc_band = bcmodel['bc_'+bd][:]

    if ((input.plx > 0.)):
        # load up bolometric correction grid
        # only K-band for now
        points = (teffgrid,logggrid,fehgrid,avgrid)
        values = bc_band
        interp = RegularGridInterpolator(points,values)

        ### Monte Carlo starts here
        
        # number of samples
        nsample=int(1e5)

        # length scale for exp decreasing vol density prior in pc
        L=1350.
        
        # maximum distance to sample (in pc)
        maxdis=1e5

        # get a rough maximum and minimum distance
        tempdis=1./input.plx
        tempdise=input.plxe/input.plx**2
        maxds=tempdis+5.*tempdise
        minds=tempdis-5.*tempdise
        
        ds=np.arange(1.,maxdis,1.)
        lh = (1./(np.sqrt(2.*np.pi)*input.plxe))*np.exp( (-1./(2.*input.plxe**2))*(input.plx-1./ds)**2)
        prior=(ds**2/(2.*L**3.))*np.exp(-ds/L)
        dis = lh*prior
        dis2=dis/np.sum(dis)
        norm=dis2/np.max(dis2)
        
        # Deal with negative and positive parallaxes differently:
        if tempdis > 0:
            # Determine maxds based on posterior:
            um=np.where((ds > tempdis) & (norm < 0.001))[0]

            # Determine minds just like maxds:
            umin=np.where((ds < tempdis) & (norm < 0.001))[0]
        else:
            # Determine maxds based on posterior, taking argmax instead of tempdis which is wrong:
            um=np.where((ds > np.argmax(norm)) & (norm < 0.001))[0]

            # Determine minds just like maxds:
            umin=np.where((ds < np.argmax(norm)) & (norm < 0.001))[0]

        if (len(um) > 0):
            maxds=np.min(ds[um])
        else:
            maxds=1e5

        if (len(umin) > 0):
            minds=np.max(ds[umin])
        else:
            minds=1.

        print 'using max distance:',maxds
        print 'using min distance:',minds
        
        ds=np.linspace(minds,maxds,nsample)
        lh = (1./(np.sqrt(2.*np.pi)*input.plxe))*np.exp( (-1./(2.*input.plxe**2))*(input.plx-(1./ds))**2)
        prior=(ds**2/(2.*L**3.))*np.exp(-ds/L)
        dis = lh*prior
        dis2=dis/np.sum(dis)

        # sample distances following the discrete distance posterior
        np.random.seed(seed=10)
        dsamp=np.random.choice(ds,p=dis2,size=nsample)
        
        # Take derived additive b value from Fulton et al. (2018) from Nishiyama et al. (2008) AH/AK = 0.063 and interpolate dustmodel dataframe to determine values of reddening.
        avs = (3.1+0.063)*np.interp(x=dsamp,xp=np.concatenate(([0.0],np.array(dustmodel.columns[2:].str[3:],dtype='float'))),fp=np.concatenate(([0.0],np.array(dustmodel.iloc[0][2:]))))
        
        # NB the next line means that useav is not actually working yet
        if (useav > -99):
            avs = np.zeros(len(dsamp))+useav
        ext=avs*extfactors['a'+bd]
        
        # assume solar metallicity if no input feh is provided
        if (input.feh == -99.):
            feh = 0.0
        else:
            feh = input.feh
    
        if (input.teff == -99.):
            if ((input.jmag > -99.) & (input.kmag > -99.)):
                jkmag = (input.jmag-np.median(avs*extfactors['aj'])) - (input.kmag-np.median(avs*extfactors['ak']))
                input.teff=casagrande_jk(jkmag,feh)
            if ((input.bmag > -99.) & (input.vmag > -99.)):
                bvmag = (input.bmag-np.median(avs*extfactors['ab'])) - (input.vmag-np.median(avs*extfactors['av']))
                input.teff=casagrande_bv(bvmag,feh)
            if ((input.btmag > -99.) & (input.vtmag > -99.)):
                bvtmag = (input.btmag-np.median(avs*extfactors['abt'])) - (input.vtmag-np.median(avs*extfactors['avt']))
                input.teff=casagrande_bvt(bvtmag,feh)
                #pdb.set_trace()
            input.teffe=100.
        #else:
        #    teff=input.teff
        #    teffe=input.teffe

        np.random.seed(seed=11)
        teffsamp=input.teff+np.random.randn(nsample)*input.teffe
        
        # hack to avoid crazy Teff samples
        teffsamp[teffsamp < 1000.]=1000.
            
        map=input.mag
        mape=input.mage
        np.random.seed(seed=12)
        map_samp=map+np.random.randn(nsample)*mape
        
        # NB no extinction correction here yet since it is either:
        #   - already taken into account in ATLAS BCs below
        #   - corrected for M dwarfs further below
        absmag = -5.*np.log10(dsamp)+map_samp+5.
        
        #pdb.set_trace()
        
        # if no logg is provided, take guess from absolute mag-logg fit to solar-metallicity MIST isochrones
        # NB these coeffs are dodgy in Mv, but pretty good in Mk
        if (input.logg == -99.):
            
            if ((band == 'vmag') | (band == 'vtmag')):
                fitv=np.poly1d([ 0.00255731, -0.07991211,  0.85140418,  1.82465197]) 
                input.logg=fitv(np.median(absmag-ext))
                print 'no input logg provided, guessing (using Mv):', input.logg
                #pdb.set_trace()
            # should really be done filter by filter with a dictionary; TODO 
            else:            
                fitk=np.poly1d([-0.01234736,  0.36684517,  3.1477089 ])
                input.logg=fitk(np.median(absmag-ext))
                print 'no input logg provided, guessing (using Mk):', input.logg 
                        
        # ATLAS BCs are inaccurate for M dwarfs; use Mann et al. 2015 Mks-R relation instead
        if ((input.teff < 4100.) & (np.median(absmag-ext) > 0.)):
            if (input.feh > -99.):
                rad = 1.9305-0.3466*(absmag-ext)+0.01647*(absmag-ext)**2*(1.+0.04458*input.feh)
            else:
                rad = 1.9515-0.3520*(absmag-ext)+0.01680*(absmag-ext)**2
            
            # add 2.5% scatter in Mks-R relation 
            rad2 = rad + np.random.randn(len(rad))*np.median(rad)*0.025
		        
            lum = rad**2 * (teffsamp/teffsun)**4
		    
        # for everything else, interpolate ATLAS BCs
        else:
            if (input.teff < np.min(teffgrid)):
                return out
            if (input.teff > np.max(teffgrid)):
                return out
            if ((input.logg > -99.) & (input.logg < np.min(logggrid))):
                return out
            if ((input.logg > -99.) & (input.logg > np.max(logggrid))):
                return out
            if ((input.feh > -99.) & (input.feh < np.min(fehgrid))):
                return out
            if ((input.feh > -99.) & (input.feh > np.max(fehgrid))):
                return out
            fix=np.where(avs > np.max(avgrid))[0]
            avs[fix]=np.max(avgrid)
            fix=np.where(avs < np.min(avgrid))[0]
            avs[fix]=np.min(avgrid)

            if ((input.teff > -99.) & (input.logg > -99.)):
                #bc = interp(np.array([input.teff,input.logg,input.feh,0.]))[0]
                arr=np.zeros((len(avs),4))
                arr[:,0]=np.zeros(len(avs))+input.teff
                arr[:,1]=np.zeros(len(avs))+input.logg
                arr[:,2]=np.zeros(len(avs))+feh
                arr[:,3]=np.zeros(len(avs))+avs
                um=np.where(arr[:,3] < 0.)[0]
                arr[um,3]=0.
                bc=interp(arr)	    

                Mvbol = absmag + bc
                lum = 10**((Mvbol-Msun)/(-2.5))
                t = teffsamp/teffsun
                rad = (lum*t**(-4.))**0.5


        #pdb.set_trace()

        '''
        out.lum=np.median(lum)
        out.lumep=np.percentile(lum,84.1)-out.lum
        out.lumem=out.lum-np.percentile(lum,15.9)

        out.rad=np.median(rad)
        out.radep=np.percentile(rad,84.1)-out.rad
        out.radem=out.rad-np.percentile(rad,15.9)

        out.dis=np.median(dsamp)
        out.disep=np.percentile(dsamp,84.1)-out.dis
        out.disem=out.dis-np.percentile(dsamp,15.9)
        
        out.avs=np.median(avs)
        out.avsep=np.percentile(avs,84.1)-out.avs
        out.avsem=out.avs-np.percentile(avs,15.9)
        '''
        
        out.rad,out.radep,out.radem=getstat(rad)
        out.mabs,out.mabsep,out.mabsem=getstat(absmag-ext)
        out.lum,out.lumep,out.lumem=getstat(lum)
        out.dis,out.disep,out.disem=getstat(dsamp)
        out.avs,out.avsep,out.avsem=getstat(avs)
        #pdb.set_trace()
        out.teff=input.teff
        out.teffe=input.teffe
        out.teffep=input.teffe
        out.teffem=input.teffe
        out.logg=input.logg
        out.logge=input.logge
        out.loggep=input.logge
        out.loggem=input.logge
        out.feh=input.feh
        out.fehe = input.fehe
        out.fehep=input.fehe
        out.fehem=input.fehe
        out.plx=input.plx
        out.plxe=input.plxe

        if (plot == 'i'): # For interactive plotting
            plt.ion()
            plt.clf()
            plt.subplot(3,2,1)
            plt.hist(teffsamp,bins=100)
            plt.title('Teff')

            plt.subplot(3,2,2)
            plt.hist(lum,bins=100)
            plt.title('Lum')

            plt.subplot(3,2,3)
            plt.hist(rad,bins=100)
            plt.title('Rad')

            plt.subplot(3,2,4)
            plt.hist(absmag,bins=100)
            plt.title('absmag')

            plt.subplot(3,2,5)
            plt.hist(dsamp,bins=100)
            plt.title('distance')

            plt.subplot(3,2,6)
            plt.hist(avs,bins=100)
            plt.title('Av')
            plt.tight_layout()
            raw_input(':')
        
        if (plot == 1): # For non-interactive plotting. Any other number for no plotting at all.
            plt.clf()
            plt.subplot(3,2,1)
            plt.hist(teffsamp,bins=100)
            plt.title('Teff')

            plt.subplot(3,2,2)
            plt.hist(lum,bins=100)
            plt.title('Lum')

            plt.subplot(3,2,3)
            plt.hist(rad,bins=100)
            plt.title('Rad')

            plt.subplot(3,2,4)
            plt.hist(absmag,bins=100)
            plt.title('absmag')

            plt.subplot(3,2,5)
            plt.hist(dsamp,bins=100)
            plt.title('distance')

            plt.subplot(3,2,6)
            plt.hist(avs,bins=100)
            plt.title('Av')


        print '   '
        print 'teff(K):',out.teff,'+/-',out.teffe
        print 'dis(pc):',out.dis,'+',out.disep,'-',out.disem
        print 'av(mag):',out.avs,'+',out.avsep,'-',out.avsem
        print 'rad(rsun):',out.rad,'+',out.radep,'-',out.radem
        print 'lum(lsun):',out.lum,'+',out.lumep,'-',out.lumem
        print 'mabs(',band,'):',out.mabs,'+',out.mabsep,'-',out.mabsem
        
        print '-----'

    ##############################################
    # case 2: input is spectroscopy + seismology #
    ##############################################
    if ((input.dnu > -99.) & (input.teff > -99.)):
        # seismic logg, density, M and R from scaling relations; this
        # is iterated, since Dnu scaling relation correction depends
        # on M
        dmass=1.
        fdnu=1.
        dnuo=input.dnu
        oldmass=1.0
        nit=0.

        while (nit < 5):

            numaxn = input.numax/numaxsun
            numaxne = input.numaxe/numaxsun
            dnun = (dnuo/fdnu)/dnusun
            dnune = input.dnue/dnusun
            teffn = input.teff/teffsun
            teffne = input.teffe/teffsun

            out.rad = (numaxn) * (dnun)**(-2.) * np.sqrt(teffn)
            out.rade = np.sqrt( (input.numaxe/input.numax)**2. + \
                4.*(input.dnue/input.dnu)**2. + \
                0.25*(input.teffe/input.teff)**2.)*out.rad

            out.mass = out.rad**3. * (dnun)**2.
            out.masse = np.sqrt( 9.*(out.rade/out.rad)**2. + \
                4.*(input.dnue/input.dnu)**2. )*out.mass

            out.rho = rho_sun * (dnun**2.)
            out.rhoe = np.sqrt( 4.*(input.dnue/input.dnu)**2. )*out.rho

            g = g_sun * numaxn * teffn**0.5
            ge = np.sqrt ( (input.numaxe/input.numax)**2. + \
                (0.5*input.teffe/input.teff)**2. ) * g

            out.logg = np.log10(g)
            out.logge = ge/(g*np.log(10.))

            # Dnu scaling relation correction from Sharma et al. 2016
            if (dnucor == 1):
                if (input.clump == 1):
                    evstate=2
                else:
                    evstate=1
                #pdb.set_trace()
                dnu,numax,fdnu=dnumodel.get_dnu_numax(evstate,input.feh,input.teff,out.mass,out.mass,out.logg,isfeh=True)
                #print out.mass,fdnu

            dmass=abs((oldmass-out.mass)/out.mass)
            oldmass=out.mass
            nit=nit+1

        print fdnu

        #pdb.set_trace()
        out.lum = out.rad**2. * teffn**4.
        out.lume = np.sqrt( (2.*out.rade/out.rad)**2. + (4.*input.teffe/input.teff)**2. )*out.lum

        print '   '
        print 'teff(K):',input.teff,'+/-',input.teffe
        print 'feh(dex):',input.feh,'+/-',input.fehe
        print 'logg(dex):',out.logg,'+/-',out.logge
        print 'rho(cgs):',out.rho,'+/-',out.rhoe
        print 'rad(rsun):',out.rad,'+/-',out.rade
        print 'mass(msun):',out.mass,'+/-',out.masse
        print 'lum(lsun):',out.lum,'+/-',out.lume
        print '-----'

        out.teff=input.teff
        out.teffep=input.teffe
        out.teffem=input.teffe
        out.feh=input.feh
        out.fehep=input.fehe
        out.fehem=input.fehe
        out.loggep=out.logge
        out.loggem=out.logge
        out.radep=out.rade
        out.radem=out.rade
        out.rhoep=out.rhoe
        out.rhoem=out.rhoe
        out.massep=out.masse
        out.massem=out.masse
        out.lumep=out.lume
        out.lumem=out.lume

        ddis=1.
        ext=0.0
        err_=0.01
        olddis=100.0

        # pick an apparent magnitude from input
        map=-99.
        if (input.vmag > -99.):
            map = input.vmag
            mape = input.vmage
            str = 'bc_v'
            avtoext=extfactors['av']

        if (input.vtmag > -99.):
            map = input.vtmag
            mape = input.vtmage
            str = 'bc_vt'
            avtoext=extfactors['avt']

        if (input.jmag > -99.):
            map = input.jmag
            mape = input.jmage
            str = 'bc_j'
            avtoext=extfactors['aj']

        if (input.kmag > -99.):
            map = input.kmag
            mape = input.kmage
            str = 'bc_k'
            avtoext=extfactors['ak']

        if (input.gamag > -99.):
            map = input.gamag
            mape = input.gamage
            str = 'bc_ga'
            avtoext=extfactors['aga']

        # if apparent mag is given, calculate distance
        if (map > -99.):
            print 'using '+str
            print 'using coords: ',input.ra,input.dec

            # iterated since BC depends on extinction
            nit=0
            while (nit < 5):

                if (nit == 0.):
                    out.avs=0.0
                else:
                    # Take derived additive b value from Fulton et al. (2018) from Nishiyama et al. (2008) AH/AK = 0.063 and interpolate dustmodel dataframe to determine values of reddening.
                    out.avs = (3.1+0.063)*np.interp(x=dsamp,xp=np.concatenate(([0.0],np.array(dustmodel.columns[2:].str[3:],dtype='float'))),fp=np.concatenate(([0.0],np.array(dustmodel.iloc[0][2:]))))[0]

                if (useav != 0.):
                    out.avs=useav
                if (out.avs < 0.):
                    out.avs = 0.0
                ext = out.avs*avtoext

                # bolometric correction interpolated from MESA

                points = (teffgrid,logggrid,fehgrid,avgrid)
                values = bc_band
                interp = RegularGridInterpolator(points, values)

                #pdb.set_trace()
                bc = interp(np.array([input.teff,out.logg,input.feh,out.avs]))[0]
                #bc = interp(np.array([input.teff,out.logg,input.feh,0.]))[0]
                
                Mvbol = -2.5*(np.log10(out.lum))+Msun
                Mvbole = np.sqrt( (-2.5/(out.lum*np.log(10.)))**2*out.lume**2)

                Mabs = Mvbol - bc
                Mabse = np.sqrt( Mvbole**2 + err_bc**2)

                ext=0. # ext already applied in BC
                logplx = (Mabs-5.-map+ext)/5.
                logplxe = np.sqrt( (Mabse/5.)**2. + (mape/5.)**2. + (err_ext/5.)**2. )

                out.plx = 10.**logplx
                out.plxe = np.log(10)*10.**logplx*logplxe

                out.dis = 1./out.plx
                out.dise = out.plxe/out.plx**2.

                ddis=abs((olddis-out.dis)/out.dis)
                #print olddis,out.dis,ddis,ext
                olddis=out.dis
		
                nit=nit+1
                #print out.dis,out.avs


            print 'Av(mag):',out.avs
            print 'plx(mas):',out.plx*1e3,'+/-',out.plxe*1e3
            print 'dis(pc):',out.dis,'+/-',out.dise

            out.disep=out.dise
            out.disem=out.dise
            
            out.mabs=Mabs
    
    return out
    
def getstat(indat):
    p16, med, p84 = np.percentile(indat,[16,50,84])
    emed1  = med - p16
    emed2  = p84 - med
    return med,emed2,emed1
    
  
def casagrande_jk(jk,feh):
    teff = 5040./(0.6393 + 0.6104*jk + 0.0920*jk**2 - 0.0330*jk*feh + 0.0291*feh + 0.0020*feh**2)
    return teff
    
def casagrande_bv(bv,feh):
    teff = 5040./(0.5665 + 0.4809*bv -0.0060*bv**2 - 0.0613*bv*feh - 0.0042*feh - 0.0055*feh**2)
    return teff

def casagrande_bvt(bvt,feh):
    teff = 5040./(0.5839 + 0.4000*bvt -0.0067*bvt**2 - 0.0282*bvt*feh - 0.00346*feh - 0.0087*feh**2)
    return teff


class obsdata():
    def __init__(self):
    
        self.ra = -99.
        self.dec = -99.
    
        self.plx = -99.
        self.plxe = -99.
    
        self.teff = -99.
        self.teffe = -99.
        self.logg = -99.
        self.logge = -99.
        self.feh = -99.
        self.fehe = -99.

        self.mag = -99.
        self.mage = -99.

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

        self.gamag = -99.
        self.gamage = -99.
        
        self.numax = -99.
        self.numaxe = -99.
        self.dnu = -99.
        self.dnue = -99.
        
        self.clump=0.
                   
    def addspec(self,value,sigma):
        self.teff = value[0]
        self.teffe = sigma[0]
        self.logg = value[1]
        self.logge = sigma[1]
        self.feh = value[2]
        self.fehe = sigma[2]

    def addmag(self,value,sigma):
        self.mag = value[0]
        self.mage = sigma[0]
               
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

    def addgaia(self,value1,value2):
        self.gamag = value1
        self.gamage = value2
        
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
        
    def addcoords(self,value1,value2):
        self.ra = value1
        self.dec = value2
        
    def addplx(self,value,sigma):
        self.plx = value
        self.plxe = sigma
        
    def addseismo(self,value,sigma):
        self.numax = value[0]
        self.numaxe = sigma[0]
        self.dnu = value[1]
        self.dnue = sigma[1]

class resdata():
    def __init__(self):
        self.teff = 0.
        self.teffe = 0.
        self.teffep = 0.
        self.teffem = 0.
        self.logg = 0.
        self.logge = 0.
        self.loggep = 0.
        self.loggem = 0.
        self.feh = 0.
        self.fehe = 0.
        self.fehep = 0.
        self.fehem = 0.
        self.rad = 0.
        self.rade = 0.
        self.radep = 0.
        self.radem = 0.
        self.mass = 0.
        self.masse = 0.
        self.massep = 0.
        self.massem = 0.
        self.rho = 0.
        self.rhoe = 0.
        self.rhoep = 0.
        self.rhoem = 0.
        self.lum = 0.
        self.lume = 0.
        self.lumep = 0.
        self.lumem = 0.
        self.avs = 0.
        self.avse = 0.
        self.avsep = 0.
        self.avsem = 0.
        self.dis = 0.
        self.dise = 0.
        self.disep = 0.
        self.disem = 0.
        self.plx = 0.
        self.plxe = 0.
        self.plxep = 0.
        self.plxem = 0.
        
        self.mabs = 0.
        self.mabse = 0.
        self.mabsep = 0.
        self.mabsem = 0.

class extinction():
    
    def __init__(self):
    
        self.ext={"ab":1.3454449, "av":1.0, "abt":1.3454449, "avt":1.0602271, "ag":1.2348743, "ar":0.88343449, "ai":0.68095687, "az":0.48308430, "aj":0.28814896, "ah":0.18152716, "ak":0.11505195, "aga":1.2348743}
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
  
    

    '''
    mass,radius=s.get_mass_radius(evstate,logz,teff,dnu,numax)
    print,mass,radius

    raw_input(':')

    # make some samples
    nsamp=1e4
    dnun=(dnu+np.random.randn(nsamp)*dnue)/dnusun
    numaxn=(numax+np.random.randn(nsamp)*numaxe)/numaxsun
    teffs=(teff+np.random.randn(nsamp)*teffe)
    teffn=teffs/5777.

    rad_sc = (numaxn) * (dnun)**(-2.) * np.sqrt(teffn)
    mass_sc = rad_sc**3. * (dnun)**2.
    rho = rho_sun * (dnun**2.)
    g = g_sun * numaxn * teffn**0.5
    logg = np.log10(g)

    #ascii.write([teffs,rad_sc,mass_sc,rho,logg], \
    #names=['teff', 'rad','mass','rho','logg'], output='epic2113_stellarsamples.txt')

    #evstate=[1.,1.]
    #logz=[-1.0,-1.0]
    #teff=[4500.,4500.]
    #mass=[1.0,1.0]
    #logg=[2.25,2.25]
    #s=asfgrid.Seism()
    #dnu,numax,fdnu=s.get_dnu_numax(evstate,logz,teff,mass,mass,logg)
    '''

