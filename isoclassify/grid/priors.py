import numpy as np
import ephem
from scipy.stats import norm

def gaussian(x, a, b, c, d):
    val = d + (a * np.exp(-(x - b)**2 / c**2))
    return val

def fehprior(feh):
    _fehprior = (0.8/0.15*np.exp(-0.5*(feh-0.016)**2./0.15**2.)
                 +0.2/0.22*np.exp(-0.5*(feh+0.15)**2./0.22**2.))
    return _fehprior

def avprior(av_in,data,i,dust,dist):
    # conversion from E(B-V) to Av from Schlafly & Finkbeiner
    fac = 2.742
    ebv_in = av_in/fac
    dm = 5.*np.log10(dist)-5.

    ra = data['ra'][i]*np.pi/180.
    dec = data['dec'][i]*np.pi/180.
    equ = ephem.Equatorial(ra, dec, epoch=ephem.J2000)
    gal = ephem.Galactic(equ)
        
    dsquared = ((gal.lon*180./np.pi - dust['lon'])**2 
                + (gal.lat*180./np.pi - dust['lat'])**2 )
    pos = np.argmin(dsquared)
    ebv = dust['vals'][pos,:,:]
    dms =np.arange(4.0,19.5,0.5)
    ebvs = np.zeros(20)
    for j in range(0,len(ebvs)):
        ebvs[j] = np.interp(dm,dms,ebv[j,:])

    med = np.median(ebvs)
    std = np.std(ebvs)
    #print gal.lon*180./np.pi,gal.lat*180./np.pi,dm,med
    #pdb.set_trace()
    return fac * gaussian(ebv_in,1.,np.median(ebvs),np.std(ebvs),0.)


def getav(data,i,dust,dist):
    # conversion from E(B-V) to Av from Schlafly & Finkbeiner
    fac = 2.742
    dm = 5.*np.log10(dist)-5.
    ra = data['ra'][i]*np.pi/180.
    dec = data['dec'][i]*np.pi/180.
    equ = ephem.Equatorial(ra, dec, epoch=ephem.J2000)
    gal = ephem.Galactic(equ)
        
    dsquared = ((gal.lon*180./np.pi - dust['lon'])**2 
                + (gal.lat*180./np.pi - dust['lat'])**2 )
    pos = np.argmin(dsquared)

    ebv = dust['vals'][pos,:,:]
    dms = np.arange(4.0,19.5,0.5)
    ebvs = np.zeros(20)
    for j in range(0,len(ebvs)):
        ebvs[j] = np.interp(dm,dms,ebv[j,:])

    med = np.median(ebvs)
    std = np.std(ebvs)
    #print gal.lon*180./np.pi,gal.lat*180./np.pi,dm,med
    #pdb.set_trace()

    #return fac*gaussian(ebv_in,1.,np.median(ebvs),np.std(ebvs),0.)
    return med,std
