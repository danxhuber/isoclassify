import numpy as np
import matplotlib.pyplot as plt
from classify_grid import *
import os, ebf
from astropy.io import ascii
import time
 

if __name__ == '__main__':

    homedir=os.path.expanduser('~/')
    model=ebf.read(homedir+'science/models/MIST/mesa.ebf')
    model['rho']=np.log10(model['rho'])
    # do this to turn off scaling relation corrections
    #model['fdnu'][:]=1.
    
    model['avs']=np.zeros(len(model['teff']))
    model['dis']=np.zeros(len(model['teff']))
    
    #dustmodel = mwdust.Combined15()m

    #data=ascii.read('../../apokasc/tgas_apokasc_obs.txt')
    #data=ascii.read('../../apokasc/tgas_dwarfs_obs.txt')
    data=ascii.read('apokasc/tgas_combined_obs.txt')

    #f = open('../../apokasc/tgas_combined_derived_grid.txt', 'w')
    #f = open('../../apokasc/tgas_combined_derived_grid_seismo_dnucor.txt', 'w')
    #f = open('../../apokasc/tgas_combined_derived_grid_plx.txt', 'w')

    f = open('temp.dat','w')

    f.write('kic teff teffep teffem logg loggep loggem feh fehep fehem rad radep radem mass massep massem lum lumep lumem rho rhoep rhoem dis disep disem av avep avem \n')

    for i in range(0,len(data['KIC'])):
        
        #if (i < 493):
        #    continue
        
        #if (data['KIC'][i] != 5456612):
        #    continue

        #if (data['logg'][i] < 3.5):
        #    continue
        
        x=obsdata()
	x.addcoords(data['ra'][i],data['dec'][i])
        
	x.addbvt([data['btmag'][i],data['vtmag'][i]],[data['sig_btmag'][i],data['sig_vtmag'][i]])
        x.addjhk([data['jmag'][i],data['hmag'][i],data['kmag'][i]],\
                 [data['sig_jmag'][i],data['sig_hmag'][i],data['sig_kmag'][i]])

        x.addgriz([data['gmag'][i],data['rmag'][i],data['imag'][i],data['zmag'][i]],\
                [data['sig_gmag'][i],data['sig_rmag'][i],data['sig_imag'][i],data['sig_zmag'][i]])
	
        x.addspec([data['teff'][i],-99.0,data['feh'][i]],[data['sig_teff'][i],0.0,data['sig_feh'][i]])

        if (data['sig_kmag'][i] > 0.5):
            continue

        # added errors are based on scatter in apokasc sample
        #dnusig=data['sig_dnu'][i]
        #if (data['sig_dnu'][i]/data['dnu'][i] < 0.005):
        #    dnusig = data['dnu'][i]*0.005
            
        #numaxsig=data['sig_numax'][i]
        #if (data['sig_numax'][i]/data['numax'][i] < 0.01):
        #    numaxsig = data['numax'][i]*0.01
        
        dnusig=np.sqrt( data['sig_dnu'][i]**2. + (data['dnu'][i]*0.005)**2.)
        numaxsig=np.sqrt( data['sig_numax'][i]**2. + (data['numax'][i]*0.01)**2.)

        if (data['numax'][i] > 0.):
            x.addseismo([data['numax'][i],data['dnu'][i]],[numaxsig,dnusig])
        else:
            x.addseismo([-99.,data['dnu'][i]],[-99.,dnusig])
        
	    
#        pdb.set_trace()
	#plxe=data['sig_plx'][i]
        #plx=data['plx'][i]
        #plx=1e3/791.88
        #plxe=data['plx'][i]*0.001
        #x.addplx(plx/1e3,plxe/1e3)

#        print x.plxe/x.plx
       	print data['KIC'][i],x.teff,x.feh
        print 'BtVt:',x.vtmag,x.vtmage,x.btmag,x.btmage
        print 'griz:',x.gmag,x.gmage,x.rmag,x.rmage,x.imag,x.image,x.zmag,x.zmage
        print 'JHK',x.jmag,x.jmage,x.hmag,x.hmage,x.kmag,x.kmage
        print '----'
        
        #if (x.plxe/x.plx < 1.0):
        t1 = time.clock()
        paras=classify(input=x,model=model,dustmodel=0,doplot=0)
        t2 = time.clock()
        print t2-t1
                
        f.write('%10i %10.1f %8.1f %8.1f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %8.3f %8.3f %8.3f\n' % (data['KIC'][i], paras.teff, paras.teffep, paras.teffem, \
        paras.logg, paras.loggep, paras.loggem, \
        paras.feh, paras.fehep, paras.fehem, \
        paras.rad, paras.radep, paras.radem, \
        paras.mass, paras.massep, paras.massem, \
        paras.lum, paras.lumep, paras.lumem, \
        paras.rho, paras.rhoep, paras.rhoem, paras.dis, paras.disep, paras.disem, paras.avs, paras.avsep, paras.avsem))     
        #f.close()
        #pdb.set_trace()
        raw_input(':')
	    
	print '---------------------------------------------'
            

    #raw_input(':')
        
    f.close()
    
    #x.addgriz([12.379,11.614,11.378,0.],[0.04,0.04,0.04,0.0])
    #x.addplx(0.00218,0.0000327)
    
    #EPIC2113
    '''
    x.addbv([13.770,12.611],[0.020,0.031])
    x.addjhk([10.694,10.177,10.035],[0.023,0.023,0.021])
    x.addgriz([13.161,12.223,12.087,0.0],[0.031,0.040,0.23,0.0])
    x.addspec([4790.,-99.0,0.38],[90.,0.0,0.08])
    x.addseismo([223.7,16.83],[5.4,0.17])
    x.addplx(0.001270435161529427,6.458185669608908e-05)
    
    paras=classify(input=x,model=model,ebv=0.0,ebve=0.0,bcerr=0.0,doplot=1)
    '''
    
    '''
    bn=np.arange(0.,1.5,0.02)

    plt.clf()
    data=ascii.read('../apokasc/tgas_match/apokasc_tgas_kepler.txt')
    plt.hist(data['sig_plx']/data['plx'],bins=bn,label='APOKASC giants')

    data=ascii.read('../apokasc/tgas_match/dwarfs_tgas_kepler.txt')
    plt.hist(data['sig_plx']/data['plx'],bins=bn,label='dwarfs')

    plt.xlabel('sig_plx/plx')
    plt.ylabel('# of stars')

    plt.legend(loc='best',numpoints=1,handletextpad=0.25,prop={'size':19},\
    handlelength=1.0)
    '''

    '''
    um=np.where(model['feh'] == 0.0)
    plt.plot(model['hmag'][um]-model['kmag'][um],model['jmag'][um]-model['hmag'][um],'.')
    plt.plot(data['hmag']-data['kmag'],data['jmag']-data['hmag'],'.')

    model=ebf.read(homedir+'science/models/DSEP/dsep.ebf')
    um=np.where(model['feh'] == 0.0)
    plt.plot(model['hmag'][um]-model['kmag'][um],model['jmag'][um]-model['hmag'][um],'.')
    '''
