# Isochrone Stellar Classification Codes

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.573372.svg)](https://doi.org/10.5281/zenodo.573372)

Python codes to perform simple stellar classifications given any set of input observables. Details are described in [Huber et al. 2017](http://adsabs.harvard.edu/abs/2017ApJ...844..102H), [Berger et al. 2020](https://ui.adsabs.harvard.edu/abs/2020arXiv200107737B/abstract) and [Berger et al. 2023](https://ui.adsabs.harvard.edu/abs/2023arXiv230111338B/abstract) - please cite these papers when using the code.

## Installation

1. Download `mwdust` (see https://github.com/jobovy/mwdust).

2. pip install the code:
    ```bash
    pip install isoclassify
    ```
    
4. If the install succeeded,

    ```bash
    isoclassify -h
    ```

    should output

    ```none
    usage: isoclassify [-h] {run,batch,multiproc,scrape-output} ...

    optional arguments:
      -h, --help            show this help message and exit

    subcommands:
      {run,batch,multiproc,scrape-output}
        run                 run isoclassify
    ```

5. **Optional:** Set an environment variable for a path to store the MESA models downloaded in step 4. Otherwise, skip this step and use the default location.

    ```bash
    export ISOCLASSIFY=/path/to/data/dir
    ```

6. Download models into the `isoclassify` data directory created upon installation or specified in step 5.

   Bolometric correction grid:
   ```bash
   https://drive.google.com/file/d/1vbgKqrghAWJrDJnRBLMlv--J8VSVUOpi/view?usp=sharing
    ```
   MESA models:
    ```bash
   https://drive.google.com/drive/folders/1GC81YxBvMF2wu4TdL0vItvdYv-zsaLI1?usp=sharing
    ```
   Parsec models:
   ```bash
   https://drive.google.com/file/d/1gURKVsL5jPxWZ08ipzALXiUnX7lsjoGe/view?usp=sharing
   ```


## Usage

### Grid Modeling:

Derives posterior distributions for stellar parameters (teff, logg, feh, rad, mass, rho, lum, age, distance, av) through direct intregration of isochrones, given any set of observables (photometry, spectroscopy, parallax, asteroseismology) and priors. Can fit for extinction or use reddening map, depending on input parameters. The Parsec grid has been amended for M dwarfs using empirical relations by [Mann et al. 2019](https://ui.adsabs.harvard.edu/abs/2019ApJ...871...63M/abstract); see [Berger et al. 2023](https://ui.adsabs.harvard.edu/abs/2023arXiv230111338B/abstract) for details. <br />


### Direct Method:

Uses bolometric corrections and extinction maps to derive stellar parameters using direct physical relations. Masses and densities are only calculated for cool stars within the empirical relations by [Mann et al. 2019](https://ui.adsabs.harvard.edu/abs/2019ApJ...871...63M/abstract) (to derive masses for any star, use grid modeling mode described above). Options are: <br />

1.  RA/DEC + Photometry + (Spectroscopy) + Parallax -> Teff, R, L, distance, Av (+ M & rho for cool stars). Uses Monte-Carlo method to implement exponentially decreasing volume density prior with a length scale of 1.35 kpc (Astraatmadja & Bailer-Jones 2016)

1. RA/DEC + Asteroseismology + Spectroscopy -> logg, rho, rad, mass, lum, distance, Av. 

Bolometric corrections are interpolated in (Teff, logg, FeH, Av) from the MIST grid (http://waps.cfa.harvard.edu/MIST/model_grids.html)


### List of Input Parameters:

#### Parameters without uncertainties: 
ra    ... Right Ascension (degrees) <br />
dec   ... Declination (degrees) <br />
band  ... Photometric band to use for absolute magnitude (see nomenclature below) <br />
dust  ... Extinction map to use (none,allsky,green19) <br />
grid  ... Stellar evolution grid to use (mesa,parsec) <br />
dmag  ... Contrast for secondary star (mag). Must be in same bandpass as "band" <br />

#### Parameters with uncertainties (e.g. plx_err): 
plx   ... Parallax (arcseconds) <br />
teff  ... Effective Temperature (K) <br />
logg  ... Surface Gravity (cgs) <br />
feh   ... Metallicity (dex) <br />
lum   ... Luminosity (solar) <br />
bmag  ... Johnson B (mag) <br />
vmag  ... Johnson V (mag) <br />
btmag ... Tycho B (mag) <br />
vtmag ... Tycho V (mag) <br />
gmag  ... Sloan g (mag) <br />
rmag  ... Sloan r (mag) <br />
imag  ... Sloan i (mag) <br />
zmag  ... Sloan z (mag) <br />
gamag ... Gaia G (mag) <br />
bpmag ... Gaia Bp (mag) <br />
rpmag ... Gaia Rp (mag) <br />
jmag  ... 2MASS J (mag) <br />
hmag  ... 2MASS H (mag) <br />
kmag  ... 2MASS Ks (mag) <br />
numax ... Frequency of maximum power (muHz) <br />
dnu   ... Asteroseismic frequency of maximum power (muHz) <br />


### Command line interface

The preferred mode to run isoclassify is through a command line interface (CLI), which allows single star processing as well as batch processing of many stars.

```bash
isoclassify run <mode> <star name> --csv <csv file> --outdir <output directory> --plot <plotmode> 
```

1. `<mode>` direct or grid
1. `<star name>` tells isoclassify which row of <csv file> to pull the input parameters
1. `<output directory>` will contain the isoclassify output files
1. `<csv file>` contains as columns parameters that are passed to isoclassify
1. `<plotmode>` tells isoclassify whether pop interactive plots `show`, save to a png `save-png`, or to not plot at all `none`.

The file `examples/example.csv` inclides a few example stars that cover various use cases, e.g.:

```bash
mkdir -p output/sol # make sure output directory exists
isoclassify run direct sol --csv examples/example.csv --outdir output/sol --plot show
```

The CLI also makes parallel processing easy.

```bash
# generate batch scripts
isoclassify batch direct examples/example.csv -o output > isoclassify.tot 

# Run with GNU parallel
parallel :::: isoclassify.tot

# Combine outputs into one CSV file
isoclassify scrape-output 'output/*/output.csv' output.csv
```

You can also run python-based multiprocessing through joblib and memory-mapping to reduce RAM overhead.

```bash
isoclassify multiproc <mode> <num processes> <input csv> <output csv> --baseoutdir <base output directory> --plot <plot mode>
```

1. `<mode>` direct or grid
1. `<num processes>` number of processes to run at a time. -1 uses all processors, including hyperthreading for machines that have it (num cores * 2). -1 should be used with caution depending upon the # of processors and RAM in your machine, as a single process can easily consume ~4GB of RAM by itself
1. `<input csv>` location and name of csv file containing all stars you want to run
1. `<output csv>` location and name of csv file, scraping all results from folders in <base output directory>
1. `<base output directory>` location of individual stellar output folders (default = ./)
1. `<plot mode>` tells isoclassify whether to save to a png `save-png` (default), save to a pdf `save-pdf`, or to not plot at all `none`. Does not use `show` as above due to issues with multiprocessing functionality.

This will run all stars within the designated csv file and produce an output csv file, so there no need to use the scrape-output function detailed above in combination with this script.

Note: isoclassify is not natively setup to run in ipython notebooks or be imported with other scripts. The example directory includes notebooks constructed by hacking  pre-defining functions that allows direct interaction with posteriors. These are provided for guidance only and are not supported.


### Default Use Case Example

In the Gaia era, a typical use-case for isoclassify is to derive the mass, radius, luminosity, density and age for a star with existing constraints on Teff and metallicity from spectroscopy, parallax from Gaia, and photometry from various surveys (e.g. Gaia and 2MASS). The recommended procedure is to first use direct-mode to determine the luminosity directly through the bolometric flux and distance:

```bash
mkdir -p output/piMen # make sure output directory exists
isoclassify run direct piMendirect --csv examples/example.csv --outdir output/piMendirect --plot show
```
In this example, the bolometric flux is estimated through 2MASS K-band (less sensitive to extinction) and the bolometric correction grid. 

Once the luminosity is known, we can provide it as an input to grid-mode, together with the spectroscopic Teff and metallicity:
```bash
isoclassify run grid piMengrid --csv examples/example.csv --outdir output/piMengrid --plot show
```
