import os
import copy
import glob
import h5py
import sys

import numpy as np
from matplotlib import pylab as plt
from matplotlib import rcParams
import pandas as pd
import astropy.units as units
from astropy.coordinates import SkyCoord
import mwdust
import pdb

from isoclassify.direct import classify as classify_direct
from isoclassify.grid import classify as classify_grid
from isoclassify.extinction import *
from isoclassify import DATADIR

CONSTRAINTS = [
    'teff','logg','feh','lum','gmag','rmag','imag','zmag','jmag','hmag','kmag',
    'gamag','bpmag','rpmag','parallax', 'bmag','vmag', 'btmag','vtmag','numax','dnu','dmag'
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

        # Set tight layout bounds to make shown figures clearer
        fig1 = plt.figure('hrd')
        fig1.set_tight_layout(True)
        fig2 = plt.figure('posteriors')
        fig2.set_tight_layout(True)

        input('[press return to continue]:')
    elif pipe.plotmode.count('save')==1:
        pipe.savefig()

    pipe.to_csv()

def runMulti(modelGridIn, **kw):
    if not os.path.exists(kw['outdir']):
        os.makedirs(kw['outdir'])
    print('running',kw['id_starname'])
    f = open(os.path.join(kw['outdir'], 'output.log'), 'w')
    sys.stdout = f

    if kw['method']=='direct':
        pipe = PipelineDirect(**kw)
        pipe.run()
    elif kw['method']=='grid':
        pipe = PipelineGrid(**kw)
        pipe.runMulti(modelGridIn)
    else:
        assert False, "method {} not supported ".format(kw['method'])

    pipe.savefig()
    pipe.to_csv()
    sys.stdout = sys.__stdout__
    f.close()

class Pipeline(object):
    def __init__(self, **kw):
        assert 'csv' in kw, "must pass csv as keyword"
        assert 'outdir' in kw, "must pass outdir as keyword"
        #assert kw.has_key('csv'), "must pass csv as keyword"
        #assert kw.has_key('outdir'), "must pass outdir as keyword"

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
        df['id_starname'] = df['id_starname'].astype(str)

        if (len(df.id_starname.drop_duplicates())!=len(df)):
            print('dropping duplicates')
            df=df.drop_duplicates(subset='id_starname')

        df.index = df.id_starname
        star = df.loc[self.id_starname]

        self.dust = star.dust
        self.grid = star.grid

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

    def addlum(self,x):
        keys = 'lum'.split()
        val = [self.const[key] for key in keys]
        err = [self.const[key+'_err'] for key in keys]
        x.addlum(val,err)

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

    def addgaia(self,x):
        keys = 'gamag bpmag rpmag'.split()
        val = [self.const[key] for key in keys]
        err = [self.const[key+'_err'] for key in keys]
        x.addgaia(val,err)

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

    def adddmag(self,x):
        x.adddmag(self.const['dmag'], self.const['dmag_err'])

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
        #pdb.set_trace()
        print("id_starname {}".format(self.id_starname))
        print("dust:", self.dust)
        print("grid:", self.grid)
        print("band:", self.const['band'])
        
        if pd.isnull(self.dust):
            self.dust=-99

        if pd.isnull(self.grid):
            self.grid=-99

        for key in CONSTRAINTS:
            if np.isnan(self.const[key]):
                self.const[key]=-99
                self.const[key+'_err']=-99
            print(key, self.const[key], self.const[key+'_err'])

        for key in COORDS:
            if np.isnan(self.const[key]):
                self.const[key]=-99
            print(key, self.const[key])

    def savefig(self):
        labels = plt.get_figlabels()
        _, ext = self.plotmode.split('-')
        for label in plt.get_figlabels():
            fn = os.path.join(self.outdir,'{}.{}'.format(label,ext))
            fig = plt.figure(label)
            fig.set_tight_layout(True)
            plt.savefig(fn)
            print("created {}".format(fn))
            plt.close()

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
        out.to_csv(self.csvfn, header=False)
        print("created {}".format(self.csvfn))



class PipelineDirect(Pipeline):
    outputcols = {
        'dir_dis': 'dis',
        'dir_avs': 'avs',
        'dir_rad': 'rad',
        'dir_lum': 'lum',
        'dir_teff': 'teff',
        'dir_mabs': 'mabs',
        'dir_mass': 'mass',
        'dir_rho': 'rho',
        'dir_fbol': 'fbol'
    }

    def run(self):
        self.print_constraints()

        fn = os.path.join(DATADIR,'bcgrid.h5')
        bcmodel = h5py.File(fn,'r', driver='core', backing_store=False)

        dustmodel,ext = query_dustmodel_coords(self.const['ra'],self.const['dec'],self.dust)

        x = classify_direct.obsdata()
        self.addspec(x)
        #self.addlum(x)
        self.addjhk(x)
        self.addbv(x)
        self.addbvt(x)
        self.addgriz(x)
        self.addgaia(x)
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
        'iso_teffsec':'teffsec',
        'iso_radsec':'radsec',
        'iso_masssec':'masssec',
        'iso_rhosec':'rhosec',
        'iso_loggsec':'loggsec',
    }
    def run(self):
        self.print_constraints()

#        model = ebf.read(os.path.join(DATADIR,'mesa.ebf'))
        fn = os.path.join(DATADIR,self.grid+'.h5')
        modfile = h5py.File(fn,'r', driver='core', backing_store=False)
        model = {'age':np.array(modfile['age']),\
        'mass':np.array(modfile['mass']),\
        'feh_init':np.array(modfile['feh']),\
        'feh':np.array(modfile['feh_act']),\
        'teff':np.array(modfile['teff']),\
        'logg':np.array(modfile['logg']),\
        'rad':np.array(modfile['rad']),\
        'lum':np.array(modfile['rad']),\
        'rho':np.array(modfile['rho']),\
        'dage':np.array(modfile['dage']),\
        'dmass':np.array(modfile['dmass']),\
        'dfeh':np.array(modfile['dfeh']),\
        'eep':np.array(modfile['eep']),\
        'bmag':np.array(modfile['bmag']),\
        'vmag':np.array(modfile['vmag']),\
        'btmag':np.array(modfile['btmag']),\
        'vtmag':np.array(modfile['vtmag']),\
        'gmag':np.array(modfile['gmag']),\
        'rmag':np.array(modfile['rmag']),\
        'imag':np.array(modfile['imag']),\
        'zmag':np.array(modfile['zmag']),\
        'jmag':np.array(modfile['jmag']),\
        'hmag':np.array(modfile['hmag']),\
        'kmag':np.array(modfile['kmag']),\
        'bpmag':np.array(modfile['bpmag']),\
        'gamag':np.array(modfile['gamag']),\
        'rpmag':np.array(modfile['rpmag']),\
        'fdnu':np.array(modfile['fdnu']),\
        'avs':np.zeros(len(np.array(modfile['gamag']))),\
        'dis':np.zeros(len(np.array(modfile['gamag'])))}

        model['rho'] = np.log10(model['rho'])
        model['lum'] = model['rad']**2*(model['teff']/5772.)**4
        # next line turns off Dnu scaling relation corrections
        # model['fdnu'][:]=1.
        model['avs']=np.zeros(len(model['teff']))
        model['dis']=np.zeros(len(model['teff']))

        dustmodel,ext = query_dustmodel_coords(self.const['ra'],self.const['dec'],self.dust)

        # Instantiate model
        x = classify_grid.obsdata()
        self.addcoords(x)
        self.addspec(x)
        self.addlum(x)
        self.addjhk(x)
        self.addgriz(x)
        self.addgaia(x)
        self.addbv(x)
        self.addbvt(x)
        self.addseismo(x)
        self.addplx(x)
        self.adddmag(x)
        self.paras = classify_grid.classify(
            input=x, model=model, dustmodel=dustmodel,ext=ext,
            plot=self.plot, useav=0, band=self.const['band']
        )

    def runMulti(self,modelGridIn):
        self.print_constraints()

        dustmodelIn,extIn = query_dustmodel_coords(self.const['ra'],self.const['dec'],self.dust)

        # Instantiate model
        x = classify_grid.obsdata()
        self.addcoords(x)
        self.addspec(x)
        self.addlum(x)
        self.addjhk(x)
        self.addgriz(x)
        self.addgaia(x)
        self.addbv(x)
        self.addbvt(x)
        self.addseismo(x)
        self.addplx(x)
        self.adddmag(x)
        self.paras = classify_grid.classify(
            input=x, model=modelGridIn, dustmodel=dustmodelIn,ext=extIn,
            plot=self.plot, useav=0, band=self.const['band']
        )

def _csv_reader(f):
    row = pd.read_csv(f,header=None, index_col=0).squeeze()
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
            print(i)
        try:
            df.append(_csv_reader(f))
        except ValueError:
            print("{} failed".format(f))


    df = pd.DataFrame(df)
    return df
