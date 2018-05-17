import sys
import ebf
import numpy as np
import scipy.interpolate
import pdb
import asfgrid
import h5py
import pandas as pd
import multiprocessing
import time
from dustmaps.bayestar import BayestarQuery
from astropy.io import ascii
from classify_direct import *

# Define the models as global to ensure they are passed to the multiprocessing functions, in addition to the seisBool and dataframe containing all data:
global dnumodelin
global bcmodelin
global dustmodelin
global seisBool
global df

# This function is what the multiprocessing loop calls. This will initialize the obsdata class and provide the necessary inputs for classify_direct to run correctly. It returns the output from classify_direct (parameters) plus the name of the star:
def compute_properties(dataFrameIndex):
    # Time how long code runs:
    t0 = time.time()
    
    # Output the position within the df:
    print "Star",dataFrameIndex,"out of",len(df)
    
    # Locate corresponding row in df:
    starProp = df.iloc[dataFrameIndex]
    
    # Define all input parameters according to the column names within your df:
    starID = starProp['KIC']
    raIn = starProp['_RA']
    decIn = starProp['_DE']
    metallicityIn = starProp['[Fe/H]i']
    metallicityInError = starProp['e_[Fe/H]i']
    loggIn = starProp['loggi']
    loggInError = starProp['e_loggi']
    teffIn = starProp['Teffi']
    teffInError = starProp['e_Teffi']
    jMagIn = starProp['kic_jmag']
    jMagInError = starProp['e_kic_jmag']
    hMagIn = starProp['kic_hmag']
    hMagInError = starProp['e_kic_hmag']
    kMagIn = starProp['kic_kmag']
    kMagInError = starProp['e_kic_kmag']
    
    # These values will be different depending on whether you are calling the asteroseismic or parallax part of the code:
    if seisBool == False:
        parallaxSysError = 0.03 # In mas, based on Lindegren et al. (2018) value
        parallaxIn = starProp['parallax'] + parallaxSysError # Make sure parallax and errors are in units of arcseconds for classify_direct
        parallaxInError = starProp['parallax_error']
        
        # Below are some checks to make sure the code does not crash by skipping stars with problematic input values:
        if metallicityIn > 0.5:
            print "Metallicity value too high. Setting value to 0.5 dex."
            metallicityIn = 0.5
        if parallaxIn == 0.0 or parallaxInError == 0.0:
            paras = resdata()
            print "****** Zero parallax encountered for",starID
            print "Teff_in:",teffIn,"+/-",teffInError
            print "logg_in:",loggIn,"+/-",loggInError
            print "[Fe/H]_in:",metallicityIn,"+/-",metallicityInError
            print "Jmag_in:",jMagIn,"+/-",jMagInError
            print "Hmag_in:",hMagIn,"+/-",hMagInError
            print "Kmag_in:",kMagIn,"+/-",kMagInError
            print "Parallax_in",parallaxIn/1000.,"+/-",parallaxInError/1000. # In units of arcseconds
            print "Skipping this star."
        elif np.isnan(parallaxIn) or np.isnan(parallaxInError):
            paras = resdata()
            print "****** NaN parallax encountered for",starID
            print "Teff_in:",teffIn,"+/-",teffInError
            print "logg_in:",loggIn,"+/-",loggInError
            print "[Fe/H]_in:",metallicityIn,"+/-",metallicityInError
            print "Jmag_in:",jMagIn,"+/-",jMagInError
            print "Hmag_in:",hMagIn,"+/-",hMagInError
            print "Kmag_in:",kMagIn,"+/-",kMagInError
            print "Parallax_in",parallaxIn/1000.,"+/-",parallaxInError/1000. # In units of arcseconds
            print "Skipping this star."
        elif (np.isnan(jMagIn) and np.isnan(hMagIn) and np.isnan(kMagIn)) or (np.isnan(jMagInError) and np.isnan(hMagInError) and np.isnan(kMagInError)):
            paras = resdata()
            print "****** NaN photometry encountered for",starID
            print "Teff_in:",teffIn,"+/-",teffInError
            print "logg_in:",loggIn,"+/-",loggInError
            print "[Fe/H]_in:",metallicityIn,"+/-",metallicityInError
            print "Jmag_in:",jMagIn,"+/-",jMagInError
            print "Hmag_in:",hMagIn,"+/-",hMagInError
            print "Kmag_in:",kMagIn,"+/-",kMagInError
            print "Parallax_in",parallaxIn/1000.,"+/-",parallaxInError/1000. # In units of arcseconds
            print "Skipping this star."
        elif (teffIn < 1) or (np.isnan(teffIn)):
            paras = resdata()
            print "****** Bad Teff encountered for",starID
            print "Teff_in:",teffIn,"+/-",teffInError
            print "logg_in:",loggIn,"+/-",loggInError
            print "[Fe/H]_in:",metallicityIn,"+/-",metallicityInError
            print "Jmag_in:",jMagIn,"+/-",jMagInError
            print "Hmag_in:",hMagIn,"+/-",hMagInError
            print "Kmag_in:",kMagIn,"+/-",kMagInError
            print "Parallax_in",parallaxIn/1000.,"+/-",parallaxInError/1000. # In units of arcseconds
            print "Skipping this star."
        else:
            # Initialize container for relevant values:
            x=obsdata()
            
            # Add necessary values for property determination:
            x.addcoords(raIn,decIn)
            x.addspec([teffIn,loggIn,metallicityIn],[teffInError,loggInError,metallicityInError])
            x.addjhk([jMagIn,hMagIn,kMagIn],[jMagInError,hMagInError,kMagInError])
            x.addplx(parallaxIn/1000.,parallaxInError/1000.)
            
            # Compute parameters:
            print "Star ID:",starID
            print "Teff_in:",teffIn,"+/-",teffInError
            print "logg_in:",loggIn,"+/-",loggInError
            print "[Fe/H]_in:",metallicityIn,"+/-",metallicityInError
            print "Jmag_in:",jMagIn,"+/-",jMagInError
            print "Hmag_in:",hMagIn,"+/-",hMagInError
            print "Kmag_in:",kMagIn,"+/-",kMagInError
            print "Parallax_in",parallaxIn/1000.,"+/-",parallaxInError/1000. # In units of arcseconds
            paras=stparas(input=x,dnumodel=0,bcmodel=bcmodelin,dustmodel=dustmodelin,useav=0,dnucor=0,plot=0)
    
    
    else:
        numaxIn = starProp['numax']
        numaxInError = starProp['e_numax']
        dnuIn = starProp['dnu']
        dnuInError = starProp['e_dnu']
    
        if metallicityIn > 0.5:
            print "Metallicity value too high. Setting value to 0.5 dex."
            metallicityIn = 0.5
        if numaxIn < 0 or dnuIn < 0 or numaxInError < 0 or dnuInError < 0 or np.isnan(numaxIn) or np.isnan(dnuIn) or np.isnan(numaxInError) or np.isnan(dnuInError):
            paras = resdata()
            print "****** Bad asteroseismic parameters encountered for",starID
            print "Teff_in:",teffIn,"+/-",teffInError
            print "logg_in:",loggIn,"+/-",loggInError
            print "[Fe/H]_in:",metallicityIn,"+/-",metallicityInError
            print "Jmag_in:",jMagIn,"+/-",jMagInError
            print "Hmag_in:",hMagIn,"+/-",hMagInError
            print "Kmag_in:",kMagIn,"+/-",kMagInError
            print "numax_in:",numaxIn,"+/-",numaxInError
            print "dnu_in:",dnuIn,"+/-",dnuInError
            print "Skipping this star."
        elif (np.isnan(jMagIn) and np.isnan(hMagIn) and np.isnan(kMagIn)) or (np.isnan(jMagInError) and np.isnan(hMagInError) and np.isnan(kMagInError)):
            paras = resdata()
            print "****** NaN photometry encountered for",starID
            print "Teff_in:",teffIn,"+/-",teffInError
            print "logg_in:",loggIn,"+/-",loggInError
            print "[Fe/H]_in:",metallicityIn,"+/-",metallicityInError
            print "Jmag_in:",jMagIn,"+/-",jMagInError
            print "Hmag_in:",hMagIn,"+/-",hMagInError
            print "Kmag_in:",kMagIn,"+/-",kMagInError
            print "numax_in:",numaxIn,"+/-",numaxInError
            print "dnu_in:",dnuIn,"+/-",dnuInError
            print "Skipping this star."
        elif (teffIn < 1) or (np.isnan(teffIn)):
            paras = resdata()
            print "****** Bad Teff encountered for",starID
            print "Teff_in:",teffIn,"+/-",teffInError
            print "logg_in:",loggIn,"+/-",loggInError
            print "[Fe/H]_in:",metallicityIn,"+/-",metallicityInError
            print "Jmag_in:",jMagIn,"+/-",jMagInError
            print "Hmag_in:",hMagIn,"+/-",hMagInError
            print "Kmag_in:",kMagIn,"+/-",kMagInError
            print "numax_in:",numaxIn,"+/-",numaxInError
            print "dnu_in:",dnuIn,"+/-",dnuInError
            print "Skipping this star."
        else:
            # Initialize container for relevant values:
            x=obsdata()
            
            # Add necessary values for property determination:
            x.addcoords(raIn,decIn)
            x.addspec([teffIn,loggIn,metallicityIn],[teffInError,loggInError,metallicityInError])
            x.addjhk([jMagIn,hMagIn,kMagIn],[jMagInError,hMagInError,kMagInError])
            x.addseismo([numaxIn,dnuIn],[numaxInError,dnuInError])
            
            # Compute parameters:
            print "Star ID:",starID
            print "Teff_in:",teffIn,"+/-",teffInError
            print "logg_in:",loggIn,"+/-",loggInError
            print "[Fe/H]_in:",metallicityIn,"+/-",metallicityInError
            print "Jmag_in:",jMagIn,"+/-",jMagInError
            print "Hmag_in:",hMagIn,"+/-",hMagInError
            print "Kmag_in:",kMagIn,"+/-",kMagInError
            print "numax_in:",numaxIn,"+/-",numaxInError
            print "dnu_in:",dnuIn,"+/-",dnuInError
            paras=stparas(input=x,dnumodel=dnumodelin,bcmodel=bcmodelin,dustmodel=dustmodelin,useav=0,dnucor=0,plot=0)
    
    # Print time elapsed:
    t1 = time.time()
    total = t1-t0
    print "Total time elapsed for star:",total,'seconds'
    
    return starID,paras

# This function runs the outer loop of the multiprocessing code, with the number of simultaneous processes equal to the numOfThreads parameter. This cannot exceed the number of (virtual) cores on your machine.
def run_multiprocessing_pipeline(numOfThreads):
    # Initialize and run multiprocessing:
    pool = multiprocessing.Pool(processes=numOfThreads)
    print "Total number of stars:",len(df.index)
    IDs,finalPars = zip(*pool.map(compute_properties,np.arange(len(df.index))))
    
    # Return arrays containing all parameters:
    return starIDs, finalPars

############################ Code Start: ###################################################

# Start timer for how long it takes to go through all data:
timeStart = time.time()

# Read in dataframe containing all necessary parameters for radius determination:
df = pd.read_csv('Kepler_GaiaDR2_Matches_2.csv') # Use whatever file you like, just make sure you change the starProp calls above

# Define boolean to determine whether we want to use the parallax or seismology inputs:
seisBool = False

# Read in bolometric correction and dust models:
bcmodelin = h5py.File('bcgrid.h5','r',driver='core',backing_store=False)
dustmodelin = BayestarQuery(version='bayestar2017')

# Decide if the asfgrid model should be read:
if seisBool:
    dnumodelin = asfgrid.Seism()

# Call multiprocessing function that returns all computed parameters:
starIDs, finalPars = run_multiprocessing_pipeline(multiprocessing.cpu_count())

# 3. Save data into df:
finalParDF = pd.DataFrame([vars(f) for f in finalPars])
finalParDF['KIC'] = pd.Series(starIDs) # Feel free to name 'KIC' whatever you like

# 4. Output df to file:
finalParDF.to_csv('Kepler_GaiaDR2_Output_Pars_2.csv',index=False) # Also name whatever you like

# End timer and print time:
timeEnd = time.time()
print "Grand total time elapsed:",timeEnd - timeStart,'seconds'

