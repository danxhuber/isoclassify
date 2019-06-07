#! /usr/bin/env python


# --------------------------------------------------------------
#       The asfgrid is a python module to compute asteroseismic 
#       parameters for a star with given stellar parameters and vice versa. 
#       Copyright (C) 2015  Sanjib Sharma, Dennis Stello

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Affero General Public License as
#     published by the Free Software Foundation, either version 3 of the
#     License, or (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Affero General Public License for more details.

#     You should have received a copy of the GNU Affero General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

# --------------------------------------------------------------

"""
A module to compute asteroseismic parameters for 
a star with given stellar parameters and vice versa. 
Author Sanjib Sharma <bugsanjib at gmail com>
Copyright (c) 2015 Sanjib Sharma, Dennis Stello
License: AGPL see <http://www.gnu.org/licenses/>.

Data files should be in current directory.

To run as a script
$./asfgrid.py --help

To use as module
::
>>> import asfgrid
>>> evstate=[1,1]
>>> logz=[-1.97,-1.98]
>>> teff=[4659.8,4903.2]
>>> dnu=[8.81,13.1]
>>> numax=[92.36,157.3]
>>> s=asfgrid.Seism()
>>> mass,radius=s.get_mass_radius(evstate,logz,teff,dnu,numax)
>>> print mass,radius
>>> logg=s.mr2logg(mass,radius)
>>> dnu,numax,fdnu=s.get_dnu_numax(evstate,logz,teff,mass,mass,logg)
>>> print dnu,numax


"""


import sys
import ebf
import numpy as np
import scipy.interpolate


__version__ = "0.0.5"


def _tocsv(filename,data,basekey=None,keylist=None,delimiter=', '):
    """
    Write a dict or npstruct to a csv file

    """
    if type(data) == dict:
        with open(filename,'w') as fp:
            if keylist==None:
                keylist=data.keys()
            if basekey == None:
                nsize=data[keylist[0]].size
            else:
                nsize=data[basekey].size        
            keylist=[key for key in keylist if data[key].size==nsize] 
            # s=str(keylist)
            # s=s[1:-1].replace("'",'')+'\n'
            s=delimiter.join([str(key) for key in keylist])+'\n'
            fp.write(s)
            for i in range(data[keylist[0]].size):
                s=', '.join([str(data[key][i]) for key in keylist])+'\n'
                fp.write(s)
    else:
        with open(filename,'w') as fp:
            if keylist==None:
                s=str(data.dtype.names)
                s=s[1:-1].replace("'",'')+'\n'
                fp.write(s)
                for temp in data:
                    s=str(temp)
                    s=s[1:-1].replace("'",'')+'\n'
                    fp.write(s)
            else:
                s=str(keylist)
                s=s[1:-1].replace("'",'')+'\n'
                fp.write(s)
                for i in range(data[keylist[0]].size):
                    s=', '.join([str(data[key][i]) for key in keylist])+'\n'
                    fp.write(s)
    print 'Written file:',filename


class _IGrid():
    def __init__(self,data1,keys):
        data=np.resize(data1,data1.size)
        data.sort(order=keys)
        self.keys=keys
        self.vnames=[temp for temp in data.dtype.names if temp not in self.keys]

        self.points=[np.unique(data[key]) for key in self.keys]
        self.values={}        
        for vname in self.vnames:
            self.values[vname]=data[vname].reshape([point.size for point in self.points])

        self.points1=tuple([data[key] for key in self.keys])
        self.values1={}        
        for vname in self.vnames:
            self.values1[vname]=data[vname]
    
    def homogenize_arrays(self,xi):
        xj=[np.asarray(t) for t in xi]
        temp=xj[0]
        for t in xj:
            temp=temp+t
        xj=[np.zeros_like(temp)+t for t in xj]
        return xj


    def get_values(self,vname,xi,fill_value='nearest'):
        fill_value1=np.nan
        if type(xi) == list:
            xi=np.array(self.homogenize_arrays(xi)).transpose()
        t1=scipy.interpolate.interpn(self.points,self.values[vname],xi,bounds_error=False,fill_value=fill_value1)
        if fill_value == 'nearest':
            ind=np.where(np.isfinite(t1)==False)[0]
            if ind.size>0:
                print 'outside interp range',ind.size,' out of ',t1.size
                if (xi.ndim == 1)&(len(self.keys)>1):
                    xi=xi.reshape([1,xi.size])
                t1[ind]=scipy.interpolate.griddata(self.points1,self.values1[vname],xi[ind],method='nearest') 
        return t1

class Seism():
    def __init__(self,datadir=None,z_solar=0.019):

# Change this to appropriate path
        if datadir is None:
            self.datadir=''
#            self.datadir='/work1/sharma/Projects/kepler/data/dnu_grid6/'
        else:
            self.datadir=datadir


# set solar reference values
#        self.radius= 6.958e8
#        self.mass=1.99e30
#        sun logg=np.log10(100.0*6.67259e-11*1.989e30/(6.958e8*6.958e8))
        self.logg_solar=4.43796037457 # cgs unit
        self.teff_solar=5777.0        # kelvin
        self.numax_solar=3090.0 # micro Hz 3090+-30
        # cannot change this
        self.dnu_solar=135.1    # micro Hz 135.1 
        self.z_solar=z_solar    # solar metallicity value


        data1=ebf.read(self.datadir+'grid_interp1.ebf','/data')
        data2=ebf.read(self.datadir+'grid_interp2.ebf','/data')
        
        self.igrid1=_IGrid(data1,['evstate','logz','mass','logg_teff'])
        self.igrid2=_IGrid(data2,['evstate','logz','mass_nu','logg_teff'])


    def logg2r(self,logg,mass):
        """
        From logg and mass compute radius

        """

        return  np.power(10.0,((self.logg_solar-logg)*0.5))*np.sqrt(mass)

    def logg2m(self,logg,radius):
        """
        From logg and radius compute mass

        """
        return np.power(10.0,logg-self.logg_solar)*radius*radius

    def logg2numax(self,logg,teff):
        """
        From logg and teff compute numax with numax_solar=3090.0 microHz

        """
        return (self.numax_solar)*np.power(10.0,logg-self.logg_solar)/(np.power(teff/self.teff_solar,0.5))

    def numax2logg(self,numax,teff):
        """
        From logg and teff compute numax with numax_solar=3090.0 microHz

        """
        return np.log10((numax/self.numax_solar)*np.sqrt(teff/self.teff_solar))+self.logg_solar


    def mr2rho(self,mass,sradius):
        """
        From mass and radius compute rho_rho_solar

        """
        return mass/np.power(sradius,3)

    def mr2logg(self,mass,radius):
        """
        From mass and radius compute logg

        """
        return self.logg_solar+np.log10(mass/(radius*radius))

    def kappa_m(self,dnu,numax):
        """
        From dnu and numax compute kappa_m

        """
        return np.power(numax/3090.0,3.0)*np.power(dnu/135.1,-4.0)

    def kappa_r(self,dnu,numax):
        """
        Not in original
        From dnu and numax compute kappa_r

        """
        return (numax/3090.0)*np.power(dnu/135.1,-2.0)

    def mass_sc(self,dnu,numax,teff):
        """
        From dnu, numax and teff compute mass according to scaling relation
        Assumes dnu_solar=135.1 microHz and numax_solar=3090.0 microHz

        """
        return np.power(numax/3090.0,3.0)*np.power(dnu/135.1,-4.0)*np.power(teff/self.teff_solar,1.5)

    def _mass_dnu(self,dnu,logg):
        """
        From dnu, logg compute mass according to scaling relation
        Assumes dnu_solar=135.1 microHz

        """
        return np.power(10.0,3*(logg-self.logg_solar))*np.power(dnu/(135.1),-4.0)
                                          

    def _quantf(self,logg,teff):
        """
        From logg and teff compute a quantity for interpolation that 
        is almost monotonic with age

        """
        return np.log10(teff)+0.5*(np.tanh((logg-4.5)/0.25)+1)*logg*0.1


    def get_dnu_numax(self,evstate,logz,teff,mini,mass,logg,fill_value='nearest',isfeh=False):
        """
        Get average seismic parameters (dnu, numax) by interpolation on 
        a grid for a given (evstate, logz, teff, dnu, numax).
        Assumption numax_solar=3090.0 microHz. 

        Args:
             evstate (array): 1) Pre RGB 2) Post RGB 
             logz or feh (array): log(Z) log metallcity or [Fe/H]=log(Z/Z_solar)
             logz    (array): log(Z) log metallcity ([Fe/H]=log(Z/Z_solar))
             teff    (array): temperature
             mini    (array): initial mass
             mass    (array): actual mass with mass loss (mass <= mini). 
             logg    (array): log(gravity)
             fill_value     : Deafault is 'nearest', to use nearest grid points 
                              in case of input values being out of grid range. 
                              Alternatively, one can use None 
        
        Returns:
             dnu     (array): Large frequency separation (micro Hz)
             numax   (array): Central frequency of max amplitude (micro Hz)
             fdnu    (array): the correction factor for Delta nu

        """
        evstate=np.asarray(evstate)
        logz=np.asarray(logz)
        if isfeh is True:
            logz=logz+np.log10(self.z_solar)
        teff=np.asarray(teff)
        mini=np.asarray(mini)
        mass=np.asarray(mass)
        logg=np.asarray(logg)

        numax=self.logg2numax(logg,teff)
        sradius=self.logg2r(logg,mass)
        logg1=self.mr2logg(mini,sradius)
        factor=self._get_fdnu(evstate,logz,teff,mini,logg1,fill_value= fill_value)
        dnu=135.1*factor*np.power(mass,0.5)/np.power(sradius,1.5)
        if (factor.size == 1)&(numax.ndim == 0):
            return dnu[0],numax,factor[0]
        else:
            return dnu,numax,factor

    def _get_fdnu(self,evstate,logz,teff,mass,logg,fill_value='nearest'):
        evstate=np.asarray(evstate)
        logz=np.asarray(logz)
        teff=np.asarray(teff)
        mass=np.asarray(mass)
        logg=np.asarray(logg)
        return self._from_mlogg('fdnu',evstate,logz,teff,mass,logg,fill_value= fill_value)



    def _from_mlogg(self,quant,evstate,logz,teff,mini,logg,fill_value='nearest'):
        """
        The driver function to perform interpolation on the grid 
        for a given (evstate, logz, teff, mini, logg) 

        Args:
             quant   (str): name of quantity for which answer is needed. 
                            For example 'fdnu', 'age', etc 
             evstate (array): 1) Pre RGB 2) Post RGB 
             logz    (array): log(Z) log metallcity ([Fe/H]=log(Z/Z_solar))
             teff    (array): temperature
             mini    (array): initial mass
             logg    (array): log(gravity)

        """
        logz=np.asarray(logz)
        logg_teff=self._quantf(logg,teff)
        return self.igrid1.get_values(quant,[evstate,logz,mini,logg_teff],fill_value= fill_value)

    def _from_freq(self,quant,evstate,logz,teff,dnu,numax,fill_value='nearest'):
        """
        The driver function to perform interpolation on the grid 
        for a given (evstate, logz, teff, dnu, numax).

        Args:
             quant   (str): name of quantity for which answer is needed. 
                            For example 'fdnu', 'age', etc 
             evstate (array): 1) Pre RGB 2) Post RGB 
             logz    (array): log(Z) log metallcity ([Fe/H]=log(Z/Z_solar))
             teff    (array): temperature
             dnu     (array): Large frequency separation (micro Hz)
             numax   (array): Central frequency of max amplitude (micro Hz)

        """
        logz=np.asarray(logz)
        logg=self.numax2logg(numax,teff)
        mass_dnu=self._mass_dnu(dnu,logg)
        logg_teff=self._quantf(logg,teff)
        return self.igrid2.get_values(quant,[evstate,logz,mass_dnu,logg_teff],fill_value= fill_value)


    def get_mass_radius(self,evstate,logz,teff,dnu,numax,fill_value='nearest',isfeh=False):
        """
        Get mass and radius of stars by interpolation on a grid 
        for a given (evstate, logz, teff, dnu, numax).
        Assumption numax_solar=3090.0 microHz. 

        Args:
             evstate     (array): 1) Pre RGB 2) Post RGB 
             logz or feh (array): log(Z) log metallcity or [Fe/H]=log(Z/Z_solar)
             teff    (array): temperature
             dnu     (array): Large frequency separation (micro Hz)
             numax   (array): Central frequency of max amplitude (micro Hz)
             fill_value     : Deafault is 'nearest', to use nearest grid points 
                              in case of input values being out of grid range. 
                              Alternatively, one can use None 
             isfeh          : A flag with default value being False. If set to 
                              True, the second argument is considered to be 
                              [Fe/H] 

        """
        evstate=np.asarray(evstate)
        logz=np.asarray(logz)
        if isfeh is True:
            logz=logz+np.log10(self.z_solar)

        teff=np.asarray(teff)
        dnu=np.asarray(dnu)
        numax=np.asarray(numax)
        mass=self._from_freq('mass',evstate,logz,teff,dnu,numax,fill_value= fill_value)
        logg=self.numax2logg(numax,teff)
        sradius=self.logg2r(logg,mass)
        if (mass.size == 1)&(evstate.ndim == 0):
            return mass[0],sradius[0]
        else:
            return mass,sradius

def _usage():
    print "NAME:"
    print "\t asfgrid 0.0.4 - computes asteroseismic freuqncies or masses"
    print "\t Copyright (c) 2015 Sanjib Sharma and Dennis Stello "
    print "\t License: AGPL see <http://www.gnu.org/licenses/>."
    print "\t Reference: Sharma et al. 2016, ApJ,822,15 \n"
    print "USAGE:"
    print "\t asfgrid inputfile \n"
    print "DESCRIPTION:"
    print "\t Outfile name is constructed from filename with suffix .out "
    print "\t Input file should be ascii as follows"
    print "\t evstate logz     teff   dnu    numax"
    print "\t 1       -1.97    4659.8 8.81   92.36"
    print "\t 1       -1.98    4903.2 13.1   157.3 \n"
    print "\t First line must contain column names"
    print "\t Column names can be in any order but need to follow names"
    print "\t given below"
    print "OPTIONS:"
    print "\t Possible input outputs are"
    print "\t 1) (evstate, logz, teff, dnu, numax)      ->(mass,radius)"
    print "\t 2) (evstate, logz, teff, mass, logg)      ->(dnu,numax,fdnu)"
    print "\t 3) (evstate, logz, teff, mass, logg, mini)->(dnu,numax,fdnu)"
    print "\t 4) (evstate, feh, teff, dnu, numax)       ->(mass,radius)"
    print "\t 5) (evstate, feh, teff, mass, logg)       ->(dnu,numax,fdnu)"
    print "\t 6) (evstate, feh, teff, mass, logg, mini) ->(dnu,numax,fdnu)"
    print "\t Third and sixth option allows for mass loss if mass<mini \n"
    print "VARIABLES:"
    print "\t evstate     (array): 1=Pre RGB tip, 2=Post RGB tip"
    print "\t logz or feh (array): Log(Z) log metallcity or [Fe/H]=log(Z/Z_solar)."
    print "\t                      If input feh, program assumes Z_solar=0.019"
    print "\t                      to convert to LogZ in the grid."
    print "\t teff        (array): Effective temperature."
    print "\t mass        (array): Actual mass; "
    print "\t                      when written as output,"
    print "\t                      mass obtained using the dnu-scaling relation"
    print "\t                      corrected with fdnu. "
    print "\t radius      (array): Radius; corresponds to the radius obtained using"
    print "\t                      the dnu-scaling relation corrected with fdnu. "
    print "\t mini        (array): Initial mass. Useful for cases with mass loss,"
    print "\t                      when actual mass is <= mini.     "
    print "\t logg        (array): Log(gravity)."
    print "\t dnu         (array): Observed large frequency separation (micro Hz);"
    print "\t                      when written as output, it corresponds to the"
    print "\t                      radial mode frequency-based dnu from the grid. "
    print "\t numax       (array): Observed frequency of max power (micro Hz);"
    print "\t                      when written as output, it corresponds to the"
    print "\t                      scaling-based numax from the grid. "
    print "\t fdnu        (array): Correction factor for Delta nu scaling relation."

if __name__  ==  '__main__':
    if len(sys.argv) == 1:
        _usage()
    elif len(sys.argv) == 2:
        if sys.argv[1] == '-help':
            _usage()
        elif sys.argv[1] == '--help':
            _usage()
        else:
            filename=sys.argv[1]
            data1=np.genfromtxt(filename,names=True)
            if data1.ndim ==0:
                data1=np.array([data1])
            data={}
            for key in data1.dtype.names:
                data[key]=data1[key]
            keylist=list(data1.dtype.names)
            s1=set(data.keys())
            status1=1
            status2=1

            if 'feh' in s1:
                data['logz']=data['feh']+np.log10(0.019)
                s1=set(data.keys())

            for key in ['evstate','logz','teff']:
                if key not in s1:
                    print 'Following columns should be present in input file' 
                    print 'evstate, logz (or feh) and teff'
                    status1=0
                    status2=0

            for key in ['mass','logg']:
                if key not in s1:
                    status1=0

            for key in ['dnu','numax']:
                if key not in s1:
                    status2=0

            if (status1+status2) !=1:
                print 'In addition to [evstate, logz, teff],' 
                print 'only one of following should be present' 
                print '[mini, mass, logg] [mass, logg] or [dnu,numax]' 
                print('Ambiguous input')
            else:
                s=Seism()
                
            if status1 == 1:
                if 'mini' in data.keys():
                    print '(evstate, logz, teff, mass, logg, mini) -> (dnu, numax,fdnu)'
                    data['dnu'],data['numax'],data['fdnu']=s.get_dnu_numax(data['evstate'],data['logz'],data['teff'],data['mini'],data['mass'],data['logg'])
                    keylist=keylist+['dnu','numax','fdnu']
                else:
                    print '(evstate, logz, teff, mass, logg) -> (dnu, numax,fdnu)'
                    data['dnu'],data['numax'],data['fdnu']=s.get_dnu_numax(data['evstate'],data['logz'],data['teff'],data['mass'],data['mass'],data['logg'])
                    keylist=keylist+['dnu','numax','fdnu']
                _tocsv(filename+'.out',data,keylist=keylist)
            elif status2 == 1:
                print '(evstate, logz, teff, dnu, numax) -> (mass,radius)'
                data['mass'],data['radius']=s.get_mass_radius(data['evstate'],data['logz'],data['teff'],data['dnu'],data['numax'])
                keylist=keylist+['mass','radius']
                _tocsv(filename+'.out',data,keylist=keylist)




