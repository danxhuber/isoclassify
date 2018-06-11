import os 
import copy
import glob
import pdb
import h5py

import numpy as np
from matplotlib import pylab as plt
import pandas as pd
import ebf
import astropy.units as units
from astropy.coordinates import SkyCoord
from dustmaps.bayestar import BayestarWebQuery

import grid.classify_grid 
import direct.classify_direct 

DATADIR = os.environ['ISOCLASSIFY']

CONSTRAINTS = [
    'teff','logg','feh','gmag','rmag','imag','zmag','jmag','hmag','kmag',
    'parallax', 'bmag','vmag', 'btmag','vtmag'
]

COORDS = ['ra','dec']

def run(**kw):
    if kw['method']=='direct':
        pipe = PipelineDirect(**kw)    
    elif kw['method']=='grid':
        pipe = PipelineGrid(**kw) 
    else:
        assert False, "method {} not supported ".format(kw['method'])

    pipe.run()
    pipe.savefig()
    pipe.to_csv()

def load_mist():
    return model

def query_dustmodel_coords(ra,dec):
    reddenMap = BayestarWebQuery(version='bayestar2017')
    sightLines = SkyCoord(ra*units.deg,dec*units.deg,frame='icrs')
    reddenContainer = reddenMap(sightLines,mode='random_sample')
    del reddenMap # To clear reddenMap from memory
    distanceSamples = np.array([0.06309573,0.07943284,0.1,0.12589255,0.15848933,0.19952627,0.25118864,0.31622776,0.3981072,0.50118726,0.6309574,0.7943282 ,1.,1.2589258,1.5848933,1.9952621,2.511887,3.1622777,3.981073,5.011873,6.3095727,7.943284,10.,12.589258,15.848933,19.952621,25.11887,31.622776,39.81073,50.11873,63.095726])*1000. # In pc, from bayestar2017 map distance samples
    
    dustModelDF = pd.DataFrame({'ra': [ra], 'dec': [dec]})
    
    for index in xrange(len(reddenContainer)):
        dustModelDF['av_'+str(round(distanceSamples[index],6))] = reddenContainer[index]

    return dustModelDF

class Pipeline(object):
    def __init__(self, **kw):
        self.id_starname = kw['id_starname']
        self.outdir = kw['outdir']        
        if kw.has_key('csv'):
            df = pd.read_csv(kw['csv'])
            assert len(df.id_starname.drop_duplicates())==len(df)
            df.index = df.id_starname
            star = df.ix[self.id_starname]

        const = {}
        for key in CONSTRAINTS:
            if key in star:
                const[key] = star[key]
                const[key+'_err'] = star[key+'_err']
            else:
                const[key] = -99
                const[key+'_err'] = 0

        for key in COORDS:
            if key in star:
                const[key] = star[key]
            else:
                const[key] = -99


        self.const = const
        self.const['ra'] = star['ra']
        self.const['dec'] = star['dec']
        self.const['band'] = star['band']
        self.pdffn = os.path.join(self.outdir,'output.pdf')
        self.csvfn = os.path.join(self.outdir,'output.csv')

    def addspec(self,x):
        keys = 'teff logg feh'.split()
        val = [self.const[key] for key in keys]
        err = [self.const[key+'_err'] for key in keys]
        x.addspec(val,err)

    def addjhk(self,x):
        keys = 'jmag hmag kmag'.split()
        val = [self.const[key] for key in keys]
        err = [self.const[key+'_err'] for key in keys]
        x.addjhk(val,err)

    def addgriz(self,x):
        keys = 'gmag rmag imag zmag'.split()
        val = [self.const[key] for key in keys]
        err = [self.const[key+'_err'] for key in keys]
        x.addgriz(val,err)
        
    def addbvt(self,x):
        keys = 'btmag vtmag'.split()
        val = [self.const[key] for key in keys]
        err = [self.const[key+'_err'] for key in keys]
        x.addbvt(val,err)
        
    def addbv(self,x):
        keys = 'bmag vmag'.split()
        val = [self.const[key] for key in keys]
        err = [self.const[key+'_err'] for key in keys]
        x.addbv(val,err)

    def addplx(self,x):
        x.addplx(self.const['parallax'], self.const['parallax_err'])
    
    def addcoords(self,x):
        x.addcoords(self.const['ra'],self.const['dec'])
    
    def addmag(self,x):
        x.addmag([self.const[self.const['band']]],[self.const[self.const['band']+'_err']])

    def addcoords(self,x):
        x.addcoords(self.const['ra'], self.const['dec'])

    def print_constraints(self):
        print "id_starname {}".format(self.id_starname)
        for key in CONSTRAINTS:
            print key, self.const[key], self.const[key+'_err']
        print "absmag constraint:", self.const['band']

        for key in COORDS:
            print key, self.const[key]
            
    def savefig(self):
        fig = plt.gcf()
        fig.set_tight_layout(True)
        plt.savefig(self.pdffn)

    def to_csv(self):
        out = {}
        out['id_starname'] = self.id_starname
        out = dict(out, **self.const)
        for outcol,incol in self.outputcols.items():
            out[outcol] = getattr(self.paras, incol)
            out[outcol+'_err1'] = getattr(self.paras, incol+'ep')
            out[outcol+'_err2'] = -getattr(self.paras, incol+'em')

        out = pd.Series(out)
        
        # Re-ordering series
        block1 = []
        block2 = []
        block3 = []
        for col in list(out.index):
            if col.count('id_starname')==1:
                block1.append(col)
                continue
            if (col.count('iso_')==1) :
                block3.append(col)
                continue

            block2.append(col)

        out = out[block1 + block2 + block3]
        out.to_csv(self.csvfn)
        print "created {}".format(self.csvfn)

class PipelineDirect(Pipeline):
    outputcols = {
        'dir_dis': 'dis',
        'dir_avs': 'avs',
        'dir_rad': 'rad',
        'dir_lum': 'lum',
    }
    
    def run(self):
        self.print_constraints()
        bcmodel = h5py.File(DATADIR+'/direct/bcgrid.h5','r',driver='core',backing_store=False)
        dustmodel = query_dustmodel_coords(self.const['ra'],self.const['dec'])

        x = direct.classify_direct.obsdata()
        self.addspec(x)
        self.addjhk(x)
        self.addbv(x)
        self.addbvt(x)
        self.addplx(x)
        self.addcoords(x)
        self.addmag(x)
        self.paras = direct.classify_direct.stparas(
            input=x,bcmodel=bcmodel,dustmodel=dustmodel,band=self.const['band'],plot=1
        )

class PipelineGrid(Pipeline):
    outputcols = {
        'iso_teff':'teff',
        'iso_logg':'logg',
        'iso_feh':'feh',
        'iso_rad':'rad',
        'iso_mass':'mass',
        'iso_age':'age',
        'iso_dis':'dis'      
    }
    def run(self):
        self.print_constraints()

        model = ebf.read(os.path.join(DATADIR,'mesa.ebf'))
        # prelims to manipulate some model variables (to be automated soon ...)
        model['rho'] = np.log10(model['rho'])
        # next line turns off Dnu scaling relation corrections
        model['fdnu'][:]=1.
        model['avs']=np.zeros(len(model['teff']))
        model['dis']=np.zeros(len(model['teff']))

        # Instantiate model
        x = grid.classify_grid.obsdata()
        self.addspec(x)
        self.addjhk(x)
        self.addgriz(x)
        self.addbv(x)
        self.addbvt(x)
        self.addplx(x)
        self.paras = grid.classify_grid.classify(
            input=x, model=model, dustmodel=0, doplot=0, useav=0
        )

def _csv_reader(f):
    row = pd.read_csv(f,header=None,squeeze=True, index_col=0)
    return row

def scrape_csv(path):
    """
    Read in isochrones csvfiles 
    Args:
        outdir (str): where to look for isochrones.csv files
    """
    fL = glob.glob(path)
    df = []

    for i, f in enumerate(fL):
        if i%100==0:
            print i
        try:
            df.append(_csv_reader(f))
        except ValueError:
            print "{} failed".format(f)


    df = pd.DataFrame(df)
    df = df.reset_index()
    return df
