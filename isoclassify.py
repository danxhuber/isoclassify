import os 
import copy
import glob

import numpy as np
from matplotlib import pylab as plt
import pandas as pd
import h5py
import ebf
import mwdust 

import grid.classify_grid 
import direct.classify_direct 

DATADIR = os.environ['ISOCLASSIFY_DATADIR']

CONSTRAINTS = [
    'teff','logg','met','gmag','rmag','imag','zmag','jmag','hmag','kmag',
    'parallax'
]

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

        self.const = const
        self.pngfn = os.path.join(self.outdir,'output.png')
        self.csvfn = os.path.join(self.outdir,'output.csv')

    def addspec(self,x):
        keys = 'teff logg met'.split()
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

    def addplx(self,x):
        x.addplx(self.const['parallax'], self.const['parallax_err'])

    def print_constraints(self):
        print "id_starname {}".format(self.id_starname)
        for key in CONSTRAINTS:
            print key, self.const[key], self.const[key+'_err']
    
    def savefig(self):
        fig = plt.gcf()
        fig.set_tight_layout(True)
        plt.savefig(self.pngfn)

    def to_csv(self):
        out = {}
        out['id_starname'] = self.id_starname
        out = dict(out, **self.const)
        for outcol,incol in self.outputcols.items():
            out[outcol] = getattr(self.paras, incol)
            out[outcol+'_err1'] = getattr(self.paras, incol+'ep')
            out[outcol+'_err2'] = -getattr(self.paras, incol+'em')

        out = pd.Series(out)
        out.to_csv(self.csvfn)

class PipelineDirect(Pipeline):
    outputcols = {
        'iso_sdis': 'dis',
        'iso_savs': 'avs',
        'iso_srad': 'rad',
        'iso_slum': 'lum',
    }

    def run(self):
        self.print_constraints()
        bcmodel = '/Users/petigura/code/isoclassify/direct/bcgrid.h5'
        dustmodel = mwdust.Combined15()

        x = direct.classify_direct.obsdata()
        self.addspec(x)
        self.addjhk(x)
        self.addgriz(x)
        self.addplx(x)
        self.paras = direct.classify_direct.stparas(
            input=x,bcmodel=bcmodel,dustmodel=dustmodel,useav=0,plot=1
        )



class PipelineGrid(Pipeline):
    outputcols = {
        'iso_steff': 'teff',
        'iso_slogg': 'logg',
        'iso_smet': 'feh',
        'iso_srad': 'rad',
        'iso_smass': 'mass',
        'iso_sage': 'age',
        'iso_sdis': 'dis'
        
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
