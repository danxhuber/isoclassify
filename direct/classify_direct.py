# code to calculate fundamental stellar parameters and distances using 
# a "direct method", i.e. adopting a fixed reddening map and bolometric 
# corrections

import numpy as np
import asfgrid
import h5py, ephem
import mwdust
from scipy.interpolate import RegularGridInterpolator
import pdb 
import pidly
import matplotlib.pyplot as plt
from astropy.stats import knuth_bin_width as knuth

def stparas(input,dnumodel=0,bcmodel=0,dustmodel=0,dnucor=0,useav=0,plot=0):

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

    # load model if they're not passed on
    if (dnumodel == 0):
        dnumodel = asfgrid.Seism()  
    if (bcmodel == 0): 
        bcmodel = h5py.File('bcgrid.h5', 'r')
    if (dustmodel == 0.):
        dustmodel = mwdust.Combined15()

    # object containing output values
    out = resdata()

    ## extinction coefficients
    extfactors=extinction()

    ########################################################################################
    # case 1: input is parallax + colors 
    ########################################################################################
    if ((input.plx > 0.)):

        # only K-band for now
	teffgrid=np.array(bcmodel['teffgrid'])
	avgrid=np.array(bcmodel['avgrid'])
        interp = RegularGridInterpolator((np.array(bcmodel['teffgrid']),\
                np.array(bcmodel['logggrid']),np.array(bcmodel['fehgrid']),\
                np.array(bcmodel['avgrid'])),np.array(bcmodel['bc_k']))

        ### Monte Carlo starts here
        
        # number of samples
        nsample=1e5

        # length scale for exp decreasing vol density prior in pc
        L=1350.

        # get a rough maximum distance
        tempdis=1./input.plx
        tempdise=input.plxe/input.plx**2
        maxds=tempdis+5.*tempdise
        
        ds=np.arange(1.,10000,1.)
        lh = np.exp( (-1./(2.*input.plxe**2))*(input.plx-1./ds)**2)
        prior=(ds**2/(2.*L**3.))*np.exp(-ds/L)
        dis = lh*prior
        dis2=dis/np.sum(dis)
        norm=dis2/np.max(dis2)
        um=np.where((ds > tempdis) & (norm < 0.001))[0]
        if (len(um) > 0):
            maxds=np.min(ds[um])
        else:
            maxds=10000.

        print 'using max distance:',maxds
        ds=np.linspace(1.,maxds,10000)
        lh = (1./(np.sqrt(2.*np.pi)*input.plxe))*\
             np.exp( (-1./(2.*input.plxe**2))*(input.plx-1./ds)**2)
        prior=(ds**2/(2.*L**3.))*np.exp(-ds/L)
        prior=np.zeros(len(lh))+1.
        dis = lh*prior
        dis2=dis/np.sum(dis)

        # sample distances following the discrete distance posterior
        np.random.seed(seed=10)
        dsamp=np.random.choice(ds,p=dis2,size=nsample)
        
        equ = ephem.Equatorial(input.ra*np.pi/180., input.dec*np.pi/180., epoch=ephem.J2000)
        gal = ephem.Galactic(equ)
        lon_deg=gal.lon*180./np.pi
        lat_deg=gal.lat*180./np.pi

        avs = 3.1*dustmodel(lon_deg,lat_deg,dsamp/1000.)
	ext=avs*extfactors.ak
	ext=0. # already in BC
    
        if (input.teff == -99.):
            teff=casagrande(jkmag,0.0)
            teffe=100.
        else:
            teff=input.teff
            teffe=input.teffe

        np.random.seed(seed=11)
        teffsamp=teff+np.random.randn(nsample)*teffe
            
        map=input.kmag
        mape=input.kmage
        np.random.seed(seed=12)
        map_samp=map+np.random.randn(nsample)*mape
        absmag = -5.*np.log10(dsamp)-ext+map_samp+5.

        if (input.teff < np.min(teffgrid)):
            return out
        if (input.teff > np.max(teffgrid)):
            return out
        #if (out.av > np.max(avgrid)):
        #    return out
        #if (out.av < np.min(avgrid)):
        #    return out

        if ((input.teff > -99.) & (input.logg > -99.) & (input.feh > -99.)):
            #bc = interp(np.array([input.teff,input.logg,input.feh,0.]))[0]
            arr=np.zeros((len(avs),4))
            arr[:,0]=np.zeros(len(avs))+input.teff
            arr[:,1]=np.zeros(len(avs))+input.logg
            arr[:,2]=np.zeros(len(avs))+input.feh
            arr[:,3]=np.zeros(len(avs))+avs
            um=np.where(arr[:,3] < 0.)[0]
            arr[um,3]=0.
            #pdb.set_trace()
            bc=interp(arr)
            #pdb.set_trace()
            
            #pdb.set_trace()
            

        Mvbol = absmag + bc	
        lum = 10**((Mvbol-Msun)/(-2.5))

        t = teffsamp/teffsun
        rad = (lum*t**(-4.))**0.5

        #pdb.set_trace()

        out.teff=input.teff
        out.teffe=input.teffe

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
        '''

        out.avs=np.median(avs)
        out.avsep=np.percentile(avs,84.1)-out.avs
        out.avsem=out.avs-np.percentile(avs,15.9)
        
        #pdb.set_trace()
        
        out.rad,out.radep,out.radem,radbn=getstat(rad)
        out.lum,out.lumep,out.lumem,lumbn=getstat(lum)
        out.dis,out.disep,out.disem,disbn=getstat(dsamp)
        #out.avs,out.avsep,out.avsem=getstat(avs)
        #pdb.set_trace()
        out.teff=input.teff
        out.teffep=input.teffe
        out.teffem=input.teffe
        out.logg=input.logg
        out.loggep=input.logge
        out.loggem=input.logge
        out.feh=input.feh
        out.fehep=input.fehe
        out.fehem=input.fehe

        if (plot == 1):
            plt.ion()
            plt.clf()
            plt.subplot(3,2,1)
            plt.hist(teffsamp,bins=100)
            plt.title('Teff')

            plt.subplot(3,2,2)
            plt.hist(lum,bins=lumbn)
            plt.title('Lum')

            plt.subplot(3,2,3)
            plt.hist(rad,bins=radbn)
            plt.title('Rad')

            plt.subplot(3,2,4)
            plt.hist(absmag,bins=100)
            plt.title('absmag')

            plt.subplot(3,2,5)
            plt.hist(dsamp,bins=disbn)
            plt.title('distance')

            plt.subplot(3,2,6)
            plt.hist(avs,bins=100)
            plt.title('Av')

        #pdb.set_trace()

        print '   '
        print 'teff(K):',out.teff,'+/-',out.teffe
        print 'dis(pc):',out.dis,'+',out.disep,'-',out.disem
        print 'av(mag):',out.avs,'+',out.avsep,'-',out.avsem
        print 'rad(rsun):',out.rad,'+',out.radep,'-',out.radem
        print 'lum(lsun):',out.lum,'+',out.lumep,'-',out.lumem
        
        print '-----'
        #raw_input(':')
        #pdb.set_trace()


    ########################################################################################
    # case 1: input is spectroscopy + seismology
    ########################################################################################
    if ((input.dnu > -99.) & (input.teff > -99.)):

        # seismic logg, density, M and R from scaling relations; this is iterated,
        # since Dnu scaling relation correction depends on M
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
            avtoext=extfactors.av

        if (input.vtmag > -99.):
            map = input.vtmag
            mape = input.vtmage
            str = 'bc_vt'
            avtoext=extfactors.avt

        if (input.jmag > -99.):
            map = input.jmag
            mape = input.jmage
            str = 'bc_j'
            avtoext=extfactors.aj

        if (input.kmag > -99.):
            map = input.kmag
            mape = input.kmage
            str = 'bc_k'
            avtoext=extfactors.ak

        if (input.gamag > -99.):
            map = input.gamag
            mape = input.gamage
            str = 'bc_ga'
            avtoext=extfactors.aga

        # if apparent mag is given, calculate distance
        if (map > -99.):
            print 'using '+str
            print 'using coords: ',input.ra,input.dec

            equ = ephem.Equatorial(input.ra*np.pi/180., input.dec*np.pi/180., epoch=ephem.J2000)
            gal = ephem.Galactic(equ)
            lon_deg=gal.lon*180./np.pi
            lat_deg=gal.lat*180./np.pi

            # iterated since BC depends on extinction
            nit=0
            while (nit < 5):

                # bolometric correction interpolated from MESA
                interp = RegularGridInterpolator((np.array(bcmodel['teffgrid']),\
                    np.array(bcmodel['logggrid']),np.array(bcmodel['fehgrid']),\
                np.array(bcmodel['avgrid'])),np.array(bcmodel[str]))

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

                out.avs = 3.1*dustmodel(lon_deg,lat_deg,out.dis/1000.)[0]

                if (useav != 0.):
                    out.avs=useav

                if (out.avs < 0.):
                    out.avs = 0.0

                #pdb.set_trace()
                ext = out.avs*avtoext
                nit=nit+1

            #pdb.set_trace()
            print 'Av(mag):',out.avs
            print 'plx(mas):',out.plx*1e3,'+/-',out.plxe*1e3
            print 'dis(pc):',out.dis,'+/-',out.dise

            out.disep=out.dise
            out.disem=out.dise
    
    return out
    
def getstat(indat):
    bn1,bn2=knuth(indat,return_bins=True)
    #(yax, xax, patches)=plt.hist(indat,bins=bn2)
    yax, xax = np.histogram(indat, bins=bn2)
    yax=yax.astype(float)
    xax=xax[0:len(xax)-1]+bn1/2.
    yax=yax/np.sum(yax)
    cprob = np.cumsum(yax)
    pos = np.argmax(yax)  
    med = xax[pos]
    temp = cprob[pos]
    ll = temp-temp*0.683
    ul = temp+(1. - temp)*0.683
    pos = np.argmin(np.abs(cprob-ll))
    emed1 = abs(med-xax[pos])
    pos = np.argmin(np.abs(cprob-ul))
    emed2 = abs(xax[pos]-med)
    #pdb.set_trace()
    
    #plt.plot([med,med],[0,np.max(yax)])

    return med,emed2,emed1,bn2
    
  
def casagrande(jk,feh):
    teff = 5040./(0.6393 + 0.6104*jk + 0.0920*jk**2 - 0.0330*jk*feh + 0.0291*feh + 0.0020*feh**2)
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

