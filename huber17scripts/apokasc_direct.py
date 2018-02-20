import sys
import ebf
import numpy as np
import scipy.interpolate
import pdb
import asfgrid
from astropy.io import ascii
from classify_direct import *
#from classify_direct import *

if __name__ == '__main__':

    dnumodel = asfgrid.Seism()  
    bcmodel = h5py.File('bcgrid.h5', 'r')
    dustmodel = mwdust.Green15()

    #data=ascii.read('../../apokasc/tgas_combined_cor_obs.txt')
    data=ascii.read('tgas_combined_obs.txt')

    #f = open('../../apokasc/tgas_combined_derived_direct_plx_lprior_useav.txt', 'w')
    #f = open('../../apokasc/tgas_combined_derived_direct_seismo.txt', 'w')
    #f = open('../../apokasc/tgas_combined_derived_direct_seismo_dnucor.txt', 'w')
    
    #f = open('../../apokasc/tgas_combined_cor_derived_direct_plx.txt', 'w')
    #f = open('../../apokasc/tgas_combined_cor_derived_direct_seismo.txt', 'w')
    #f = open('../../apokasc/tgas_combined_cor_derived_direct_seismo_dnucor.txt', 'w')
    
    

    #data=ascii.read('../../apokasc/tgas_dwarfs_spc_obs.txt')
    #f = open('../../apokasc/tgas_dwarfs_apo_derived_direct_seismo_dnucor.txt', 'w')
    #f = open('../../apokasc/tgas_dwarfs_spc_derived_direct_plx_lprior.txt', 'w')
    
    #f = open('../../apokasc/tgas_giants_payne_derived_direct_seismo.txt', 'w')
    #f = open('../../apokasc/tgas_dwarfs_spc_derived_direct_plx_lprior.txt', 'w')

    #data=ascii.read('../../apokasc/tgas_combined_obs.txt')
    #f = open('../../apokasc/tgas_combined_derived_direct_seismo3.txt', 'w')

    # read Av from grid-modeling here
    # av=ascii.read('../../apokasc/tgas_combined_derived_grid_seismo.txt')


    #f = open('../../mabsk.txt', 'w')
    f = open('../../temp.txt', 'w')

    f.write('kic teff teffep teffem logg loggep loggem feh fehep fehem rad radep radem mass massep massem lum lumep lumem rho rhoep rhoem dis disep disem av avep avem \n')
    #f.write('kic teff logg mass dis mabs \n')
    
    for i in range(0,len(data['KIC'])):

        print '---------------------------------'
        print 'KIC ',data['KIC'][i]
        print data['teff'][i],data['logg'][i],data['feh'][i]
        print data['plx'][i]/1e3,data['sig_plx'][i]/1e3
        print data['sig_plx'][i]/data['plx'][i]
        print data['jmag'][i],data['hmag'][i],data['kmag'][i]
        print data['sig_jmag'][i],data['sig_hmag'][i],data['sig_kmag'][i]

        if (data['sig_kmag'][i] > 0.5):
            #pdb.set_trace()
            continue
        
        #if (data['KIC'][i] != 5033245.):
        #    continue
        #if (data['logg'][i] < 3.5):
        #    continue
        
        # use Av from seismo
        #um=np.where(data['KIC'][i] == av['kic'])[0]
        #thisav=av['av'][um]
        #pdb.set_trace()

        x=obsdata()
        #x.addbvt([data['btmag'][i],data['vtmag'][i]],[data['sig_btmag'][i],data['sig_vtmag'][i]])

        x.addcoords(data['ra'][i],data['dec'][i])

        x.addjhk([data['jmag'][i],data['hmag'][i],data['kmag'][i]],\
                 [data['sig_jmag'][i],data['sig_hmag'][i],data['sig_kmag'][i]])


        # uncomment this for Gaia -> seismo
        '''
        x.addspec([data['teff'][i],data['logg'][i],data['feh'][i]],\
                  [data['sig_teff'][i],data['sig_logg'][i],data['sig_feh'][i]])
                  
        x.addplx((data['plx'][i])/1e3,data['sig_plx'][i]/1e3)
        if (data['plx'][i] < 0.):
            continue
        '''

        # uncomment this for seismo -> Gaia
        #'''
        x.addspec([data['teff'][i],-99.0,data['feh'][i]],[data['sig_teff'][i],0.0,data['sig_feh'][i]])

        if (data['numax'][i] == 0.):
            continue
        if (data['teff'][i] == 0.):
            continue   

        #dnusig=data['sig_dnu'][i]
        #if (data['sig_dnu'][i]/data['dnu'][i] < 0.005):
        #    dnusig = data['dnu'][i]*0.005
        dnusig=np.sqrt( data['sig_dnu'][i]**2. + (data['dnu'][i]*0.005)**2.)
            
        #numaxsig=data['sig_numax'][i]
        #if (data['sig_numax'][i]/data['numax'][i] < 0.01):
        #    numaxsig = data['numax'][i]*0.01
        numaxsig=np.sqrt( data['sig_numax'][i]**2. + (data['numax'][i]*0.01)**2.)

        if (data['numax'][i] > 0.):
            x.addseismo([data['numax'][i],data['dnu'][i]],[numaxsig,dnusig])
        else:
            x.addseismo([-99.,data['dnu'][i]],[-99.,dnusig])

        if (data['clump'][i] == 1):
            x.clump=1
        else:
            x.clump=0
        #'''
        
        paras=stparas(input=x,dnumodel=dnumodel,bcmodel=bcmodel,dustmodel=dustmodel,\
                      useav=0.,dnucor=0,plot=0)

        
        f.write('%10i %10.1f %8.1f %8.1f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %8.3f %8.3f %8.3f\n' % (data['KIC'][i], paras.teff, paras.teffep, paras.teffem, \
        paras.logg, paras.loggep, paras.loggem, \
        paras.feh, paras.fehep, paras.fehem, \
        paras.rad, paras.radep, paras.radem, \
        paras.mass, paras.massep, paras.massem, \
        paras.lum, paras.lumep, paras.lumem, \
        paras.rho, paras.rhoep, paras.rhoem, paras.dis, paras.disep, paras.disem, paras.avs, paras.avsep, paras.avsem)) 
        

        #f.write('%10i %10.1f %8.3f %8.3f %8.3f %8.3f \n' % (data['KIC'][i], paras.teff, paras.logg, paras.mass, paras.dis, paras.mabs)) 
        
           
        #raw_input(':')

    f.close()

    '''
    x.addbvt([9.,8.],[0.01,0.01])
    x.addspec([4790.,-99.,0.38],[100.,-99.,0.1])
    x.addseismo([223.7,16.8],[5.5,0.17])
    #pdb.set_trace()
    paras=stparas(input=x,dnumodel=dnumodel,bcmodel=bcmodel,dustmodel=dustmodel)
    pdb.set_trace()
    '''
