import os 
import copy
import glob
import h5py

import numpy as np
from matplotlib import pylab as plt
import pandas as pd
import ebf

from isoclassify.direct import classify as classify_direct
from isoclassify.grid import classify as classify_grid
from isoclassify.extinction import *
from isoclassify import DATADIR

CONSTRAINTS = [
    'teff','logg','feh','gmag','rmag','imag','zmag','jmag','hmag','kmag',
    'parallax', 'bmag','vmag', 'btmag','vtmag','numax','dnu'
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
    if pipe.plotmode=='show':
        plt.ion()
        plt.show()
        raw_input('[press return to continue]:')
    elif pipe.plotmode.count('save')==1:
        pipe.savefig()

    pipe.to_csv()

class Pipeline(object):
    def __init__(self, **kw):
        assert kw.has_key('csv'), "must pass csv as keyword"
        assert kw.has_key('outdir'), "must pass outdir as keyword"

        self.plotmode = kw['plot']

        # create plot (both interactive and save modes)
        if self.plotmode=='none':
            self.plot = 0 
        else:
            self.plot = 1

        self.id_starname = kw['id_starname']
        self.outdir = kw['outdir']

        # Read in inputs
        df = pd.read_csv(kw['csv'])
        assert len(df.id_starname.drop_duplicates())==len(df)
        df.index = df.id_starname
        star = df.ix[self.id_starname]

        self.dust = star.dust

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
        
    def addseismo(self,x):
        keys = 'numax dnu'.split()
        val = [self.const[key] for key in keys]
        err = [self.const[key+'_err'] for key in keys]
        x.addseismo(val,err)
    
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
        x.addmag(
            [self.const[self.const['band']]],
            [self.const[self.const['band']+'_err']]
        )

    def addcoords(self,x):
        x.addcoords(self.const['ra'], self.const['dec'])

    def print_constraints(self):
        print "id_starname {}".format(self.id_starname)
        print "dust:", self.dust
        for key in CONSTRAINTS:
            print key, self.const[key], self.const[key+'_err']

        for key in COORDS:
            print key, self.const[key]
            
    def savefig(self):
        labels = plt.get_figlabels()
        _, ext = self.plotmode.split('-')
        for label in plt.get_figlabels():
            fn = os.path.join(self.outdir,'{}.{}'.format(label,ext))
            fig = plt.figure(label)
            fig.set_tight_layout(True)
            plt.savefig(fn)
            print "created {}".format(fn)

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
        'dir_teff': 'teff',
        'dir_mabs': 'mabs',
    }
    
    def run(self):
        self.print_constraints()

        fn = os.path.join(DATADIR,'bcgrid.h5')
        bcmodel = h5py.File(fn,'r', driver='core', backing_store=False)
        
        if self.dust == 'allsky':
            dustmodel = query_dustmodel_coords_allsky(
                self.const['ra'],self.const['dec']
            )
            ext = extinction('cardelli')

        if self.dust == 'green18':
            dustmodel = query_dustmodel_coords(
                self.const['ra'],self.const['dec']
            )
            ext = extinction('schlafly16')

        if self.dust == 'none':
            dustmodel = 0
            ext = extinction('cardelli')

        x = classify_direct.obsdata()
        self.addspec(x)
        self.addjhk(x)
        self.addbv(x)
        self.addbvt(x)
        self.addplx(x)
        self.addcoords(x)
        self.addmag(x)
        self.paras = classify_direct.stparas(
            input=x, bcmodel=bcmodel, dustmodel=dustmodel, 
            band=self.const['band'], ext=ext, plot=1
        )

class PipelineGrid(Pipeline):
    outputcols = {
        'iso_age':'age',
        'iso_avs':'avs',
        'iso_dis':'dis',
        'iso_feh':'feh',
        'iso_mass':'mass',
        'iso_rad':'rad',
        'iso_lum':'lum',
        'iso_logg':'logg',
        'iso_rho': 'rho',
        'iso_teff':'teff',
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

        if self.dust == 'allsky':
            dustmodel = query_dustmodel_coords_allsky(self.const['ra'],self.const['dec'])
            ext = extinction('cardelli')
        if self.dust == 'green18':
            dustmodel = query_dustmodel_coords(self.const['ra'],self.const['dec'])
            ext = extinction('schlafly16')
        if self.dust == 'none':
            dustmodel = 0
            ext = extinction('cardelli')
            
        # Instantiate model
        x = classify_grid.obsdata()
        self.addcoords(x)
        self.addspec(x)
        self.addjhk(x)
        self.addgriz(x)
        self.addbv(x)
        self.addbvt(x)
        self.addseismo(x)
        self.addplx(x)
        self.paras = classify_grid.classify(
            input=x, model=model, dustmodel=dustmodel,ext=ext, 
            plot=self.plot, useav=0
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
