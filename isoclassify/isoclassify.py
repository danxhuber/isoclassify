import os
import glob
from argparse import ArgumentParser
from collections import OrderedDict

import pylab as pl
import numpy as np
import pandas as pd

import isoclassify.pipeline
import isoclassify.extinction
from isoclassify import DATADIR

import h5py
import pdb
from joblib import Parallel, delayed, load, dump
import gc

def main():
    p = ArgumentParser()
    subp = p.add_subparsers(title="subcommands", dest='subcommand')

    p2 = subp.add_parser('run',help='run isoclassify')
    p2.add_argument('method',choices=['direct', 'grid'],help='method')
    p2.add_argument('id_starname', help='name of star')
    p2.add_argument('-o','--outdir',default='./')
    p2.add_argument('--csv', help='path to csv file')
    p2.add_argument(
        '--plot', choices=['none', 'show', 'save-pdf','save-png'],
        default='save-png',
        help="plotmode: none: Don't create plots, "
        "show: create plot in interactive mode, "
        "save-pdf: save plot as pdf file"
    )
    p2.set_defaults(func=run)

    p2 = subp.add_parser('batch', description="Create batch jobs")
    p2.add_argument('method',choices=['direct', 'grid'],help='method')
    p2.add_argument('csv', type=str, help='list of parameters')
    p2.add_argument(
        '-o','--baseoutdir',default='./',help='base path of output directory'
    )
    p2.add_argument('--plot', default='save-png',help="passed to run command")
    p2.set_defaults(func=batch)

    p2 = subp.add_parser('multiproc', description="Run python-based multiprocessing.")
    p2.add_argument('method',choices=['direct', 'grid'],help='method')
    p2.add_argument('cpuCount', type=int, help='number of cpu cores to use for processing')
    p2.add_argument('csv', type=str, help='list of parameters')
    p2.add_argument('csvfn', help='Name of combined csv file')
    p2.add_argument(
        '-o','--baseoutdir',default='./',help='base path of output directory'
    )
    p2.add_argument(
        '--plot', choices=['none', 'save-pdf','save-png'],
        default='save-png',
        help="plotmode: none: Don't create plots, "
        "save-pdf: save plot as pdf file"
    )
    p2.set_defaults(func=run_multi)

    p2 = subp.add_parser('scrape-output', description="Combine batch outputs")
    p2.add_argument('baseoutdir', type=str,help='base path of output directory')
    p2.add_argument('csvfn', help='Name of combined csv file')
    p2.set_defaults(func=scrape_output)
    args = p.parse_args()
    args.func(args)

def run(args):
    isoclassify.pipeline.run(**vars(args))

def batch(args):
    df = pd.read_csv(args.csv,engine='python')
    print("#!/bin/sh")
    for i, row in df.iterrows():
        fmt = dict(**row)
        fmt = dict(fmt,**vars(args))
        fmt['outdir'] = "{baseoutdir:}/{id_starname:}".format(**fmt)
        s = ""
        s+="mkdir -p {outdir:}; "
        s+="isoclassify run {method:} {id_starname:} "
        s+="--outdir {outdir:} "
        s+="--csv {csv:} "
        s+="--plot {plot:} "
        s+="&> {outdir:}/output.log"
        s = s.format(**fmt)
        print(s)

def scrape_output(args):
    path = args.baseoutdir
    df = isoclassify.pipeline.scrape_csv(path)
    df.to_csv(args.csvfn,index=False)

'''
def read_modelGrid():
    fn = os.path.join(DATADIR,'parsec.h5')
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
    model['avs']=np.zeros(len(model['teff']))
    model['dis']=np.zeros(len(model['teff']))

    return model
'''

def memmap_grid():
    # Now define tempfiles to load into memory for multiple child processes:
    temp_folder = os.path.join(DATADIR,'tmp')
    if (os.path.exists(temp_folder) == False):
        os.mkdir(temp_folder)
    filename_model = os.path.join(temp_folder, 'isoclassify_modelGrid.mmap')

    if os.path.exists(filename_model):
        # If the memmapped file exists, do not read in model grid:
        print('Memory-mapped file already exists here:',filename_model)
        print('Loading file now.')

    else:
        # Otherwise, read model grid file:
        modelGridIn = read_modelGrid()

        # Dump memory contents to memory-mapped file:
        _model = dump(modelGridIn, filename_model)

        # Delete memory contents with garbage collector:
        _model = gc.collect()

    # Now load memory mapped files:
    modelGridIn = load(filename_model, mmap_mode='r')

    return modelGridIn

def scrape_csv_multi(i,f,dfOut):
    if i%100==0:
        print(i)
    try:
        dfOut.append(pd.read_csv(f,header=None,index_col=0).squeeze())
    except ValueError:
        print("{} failed".format(f))

def run_multi(args):
    df = pd.read_csv(args.csv)
    df['outdir'] = args.baseoutdir + '/' + df['id_starname']

    if args.method == 'grid':
        modelGridIn = memmap_grid()
        # Run pipeline in parallel with memmapped model grid:
        Parallel(n_jobs=args.cpuCount, max_nbytes=None, verbose=10)(delayed(isoclassify.pipeline.runMulti)(modelGridIn, **dict(**row,**vars(args))) for i,row in df.iterrows())
    else:
        # Run pipeline in parallel. No need to use memmapped model grid because the bcmodel is small:
        modelGridIn = {}
        Parallel(n_jobs=args.cpuCount, max_nbytes=None, verbose=10)(delayed(isoclassify.pipeline.runMulti)(modelGridIn, **dict(**row,**vars(args))) for i,row in df.iterrows())

    # Now output to file using shared memory DataFrame, dfOut, for eventual output to csv file:
    print("Now scraping all output.csv files into ",args.csvfn)
    path = args.baseoutdir + '/*/output.csv'
    fL = glob.glob(path)
    dfOut = []
    Parallel(n_jobs=args.cpuCount, max_nbytes=None, require='sharedmem')(delayed(scrape_csv_multi)(i,f,dfOut) for i,f in enumerate(fL))
    dfOut = pd.DataFrame(dfOut)
    dfOut.to_csv(args.csvfn,index=False)

if __name__=="__main__":
    main()
