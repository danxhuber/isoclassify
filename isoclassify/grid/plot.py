import numpy as np
import matplotlib.pyplot as plt

import random
from .priors import *
import fnmatch
import pdb

def plotinit():
    fig1 = plt.figure('posteriors',figsize=(8,12))
    fig2 = plt.figure('hrd',figsize=(12,12))
    plt.figure('posteriors')

def plotposterior(x,y,res,err1,err2,names,j,ix,iy):
    plt.rcParams['font.size']=8
    #fig = plt.figure(figsize=(8,12))
    plt.subplot(len(names),2,ix)
    plt.plot(x,np.cumsum(y))
    plt.plot([res,res],[0,1],'r')
    plt.plot([res+err1,res+err1],[0,1],'--r')
    plt.plot([res-err2,res-err2],[0,1],'--r')
    plt.ylim([0,1])
    plt.title(names[j])
    if fnmatch.fnmatch(names[j],'*rho*'):
        plt.xscale('log')
    if fnmatch.fnmatch(names[j],'*lum*'):
        plt.xscale('log')


    plt.subplot(len(names),2,iy)
    plt.plot(x,y)
    plt.plot([res,res],[0,1],'r')
    plt.plot([res+err1,res+err1],[0,1],'--r')
    plt.plot([res-err2,res-err2],[0,1],'--r')
    plt.ylim([0,np.max(y)+np.max(y)*0.1])
    plt.title(names[j])
    if fnmatch.fnmatch(names[j],'*rho*'):
        plt.xscale('log')
    if fnmatch.fnmatch(names[j],'*lum*'):
        plt.xscale('log')

    if fnmatch.fnmatch(names[j],'*feh*'):
        xt=np.arange(-2.,1.,0.01)
        yt=fehprior(xt)
        plt.plot(xt,yt*np.max(y)/np.max(yt),'--g')

def plotposterior_sec(x,y,res,err1,err2,names,j,ix,iy):
    fig3 = plt.figure('secondary',figsize=(8,12))
    plt.figure('secondary')
    plt.subplot(len(names),2,ix)
    plt.plot(x,np.cumsum(y))
    plt.plot([res,res],[0,1],'r')
    plt.plot([res+err1,res+err1],[0,1],'--r')
    plt.plot([res-err2,res-err2],[0,1],'--r')
    plt.ylim([0,1])
    plt.title(names[j])
    if fnmatch.fnmatch(names[j],'*rho*'):
        plt.xscale('log')
    if fnmatch.fnmatch(names[j],'*lum*'):
        plt.xscale('log')


    plt.subplot(len(names),2,iy)
    plt.plot(x,y)
    plt.plot([res,res],[0,1],'r')
    plt.plot([res+err1,res+err1],[0,1],'--r')
    plt.plot([res-err2,res-err2],[0,1],'--r')
    plt.ylim([0,np.max(y)+np.max(y)*0.1])
    plt.title(names[j])
    if fnmatch.fnmatch(names[j],'*rho*'):
        plt.xscale('log')
    if fnmatch.fnmatch(names[j],'*lum*'):
        plt.xscale('log')

    if fnmatch.fnmatch(names[j],'*feh*'):
        xt=np.arange(-2.,1.,0.01)
        yt=fehprior(xt)
        plt.plot(xt,yt*np.max(y)/np.max(yt),'--g')

def plotcc_auto(model,modelSel,input,ran,umran,d,g,mag1,mag2,mag3):
    plt.plot(model[mag1][ran[d]]-model[mag2][ran[d]],\
             model[mag2][ran[d]]-model[mag3][ran[d]],\
    '.',color='blue',markersize=1,zorder=-32)
    plt.plot(model[mag1][ran[g]]-model[mag2][ran[g]],\
             model[mag2][ran[g]]-model[mag3][ran[g]],\
    '.',color='red',markersize=1,zorder=-32)
    plt.plot(modelSel[mag1][umran]-modelSel[mag2][umran],\
             modelSel[mag2][umran]-modelSel[mag3][umran],\
    '.',color='black',markersize=1,zorder=-32)

    if ((vars(input)[mag1] > -99) & (vars(input)[mag2] > -99) & (vars(input)[mag3] > -99)):
        plt.errorbar([vars(input)[mag1]-vars(input)[mag2]], [vars(input)[mag2]-vars(input)[mag3]], \
                 xerr=np.sqrt(vars(input)[mag1+'e']**2+vars(input)[mag2+'e']**2), \
                 yerr=np.sqrt(vars(input)[mag2+'e']**2+vars(input)[mag3+'e']**2),color='green',elinewidth=5)

    plt.xlabel(mag1+'-'+mag2)
    plt.ylabel(mag2+'-'+mag3)
    plt.minorticks_on()
    plt.autoscale()

def plotcm_auto(model,modelSel,input,mabs,mabse,ran,umran,d,g,mag1,mag2,absmag):
    plt.plot(model[mag1][ran[d]]-model[mag2][ran[d]],\
             model[absmag][ran[d]],'.',color='blue',markersize=1,zorder=-32)
    plt.plot(model[mag1][ran[g]]-model[mag2][ran[g]], \
             model[absmag][ran[g]],'.',color='red',markersize=1,zorder=-32)
    plt.plot(modelSel[mag1][umran]-modelSel[mag2][umran], \
             modelSel[absmag][umran],'.',color='black',markersize=1,zorder=-32)

    if ((vars(input)[mag1] > -99) & (vars(input)[mag2] > -99) & (vars(input)[absmag] > -99)):
        col = vars(input)[mag1] - vars(input)[mag2]
        cole = np.sqrt(vars(input)[mag1+'e']**2 + vars(input)[mag2+'e']**2)
        plt.errorbar([col], [mabs], xerr=cole, yerr=mabse,color='green',elinewidth=5)

    plt.xlabel(mag1+'-'+mag2)
    plt.ylabel(absmag)
    plt.gca().invert_yaxis()
    plt.minorticks_on()
    plt.autoscale()

def plothrd_auto(model,modelSel,input,mabs,mabse,ran,umran,d,g):
    if (input.numax == -99):
        plt.plot(model['teff'][ran[d]],model['logg'][ran[d]],\
                '.',color='blue',markersize=1,zorder=-32)
        plt.plot(model['teff'][ran[g]],model['logg'][ran[g]],\
                '.',color='red',markersize=1,zorder=-32)
        plt.plot(modelSel['teff'][umran],modelSel['logg'][umran],\
                '.',color='black',markersize=1,zorder=-32)
        if (vars(input)['teff'] > -99 & (vars(input)['logg'] > -99)):
            plt.errorbar([input.teff], [input.logg], xerr=np.abs(input.teffe), yerr=np.abs(input.logge),\
                 color='green',elinewidth=5)

        plt.xlabel('teff')
        plt.ylabel('logg')
        plt.xlim([10000,2800])
        plt.ylim([5.5,0])

    else:
        mod_numax=3090*(10**model['logg']/27420.)*(model['teff']/5772.)**(-0.5)
        plt.semilogy(model['teff'][ran[d]],mod_numax[ran[d]],\
                 '.',color='blue',markersize=1,zorder=-32)
        plt.plot(model['teff'][ran[g]],mod_numax[ran[g]],\
                 '.',color='red',markersize=1,zorder=-32)
        plt.plot(modelSel['teff'][umran],mod_numax[umran],\
                 '.',color='black',markersize=1,zorder=-32)
        if (vars(input)['teff'] > -99 & (vars(input)['numax'] > -99)):
            plt.errorbar([input.teff], [input.numax], xerr=np.abs(input.teffe), yerr=np.abs(input.numaxe), \
                 color='green',elinewidth=5)

        plt.xlim([10000,2800])
        plt.ylim([100000,0.1])

    plt.minorticks_on()

def plothrd(model,modelSel,um,input,mabs,mabse,ix,iy):
    plt.subplots_adjust(
        left=0.08, bottom=0.05, right=0.96, top=0.96, wspace=0.31, hspace=0.26
    )

    plt.figure('hrd')
    frac = 0.01

    # Select a fractional subset of models from the entire model array according to frac parameter:
    ran=np.array(random.sample(range(len(model['teff'])),\
    int(len(model['teff'])*frac)))
    umran = np.array(random.sample(list(um),int(max(len(um)*frac,1))))

    # Choose dwarf (d) and giant (g) models:
    d=np.where(model['logg'][ran] > 3.5)[0]
    g=np.where(model['logg'][ran] < 3.5)[0]

    # Plot color-color diagrams on top, and change bands by editing the mag strings:
    plt.subplot(2,4,1)
    plotcc_auto(model=model,modelSel=modelSel,input=input,ran=ran,umran=umran,d=d,g=g,mag1='gmag',mag2='rmag',mag3='imag')

    plt.subplot(2,4,2)
    plotcc_auto(model=model,modelSel=modelSel,input=input,ran=ran,umran=umran,d=d,g=g,mag1='rmag',mag2='imag',mag3='zmag')

    plt.subplot(2,4,3)
    plotcc_auto(model=model,modelSel=modelSel,input=input,ran=ran,umran=umran,d=d,g=g,mag1='bpmag',mag2='gamag',mag3='rpmag')

    plt.subplot(2,4,4)
    plotcc_auto(model=model,modelSel=modelSel,input=input,ran=ran,umran=umran,d=d,g=g,mag1='jmag',mag2='hmag',mag3='kmag')

    # Plot color-magnitude diagram on bottom row; change mags by editing the mag strings:
    plt.subplot(2,4,5)
    plotcm_auto(model=model,modelSel=modelSel,input=input,mabs=mabs,mabse=mabse,ran=ran,umran=umran,d=d,g=g,mag1='gmag',mag2='kmag',absmag='kmag')

    plt.subplot(2,4,6)
    plotcm_auto(model=model,modelSel=modelSel,input=input,mabs=mabs,mabse=mabse,ran=ran,umran=umran,d=d,g=g,mag1='rmag',mag2='kmag',absmag='kmag')

    plt.subplot(2,4,7)
    plotcm_auto(model=model,modelSel=modelSel,input=input,mabs=mabs,mabse=mabse,ran=ran,umran=umran,d=d,g=g,mag1='bpmag',mag2='rpmag',absmag='gamag')

    # Plot logg-teff diagram in last plot unless numax is constrained, then use it instead of logg:
    plt.subplot(2,4,8)
    plothrd_auto(model=model,modelSel=modelSel,input=input,mabs=mabs,mabse=mabse,ran=ran,umran=umran,d=d,g=g)
