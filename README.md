# Isochrone Stellar Classification Codes

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.573372.svg)](https://doi.org/10.5281/zenodo.573372)

Python codes to perform stellar classifications given any set of input observables. Details are described in [Huber et al. 2017](http://adsabs.harvard.edu/abs/2017ApJ...844..102H) and [Berger et al. 2020](https://ui.adsabs.harvard.edu/abs/2020arXiv200107737B/abstract) - please cite these papers when using the code.

## Installation

```bash
# Clone repo
git clone https://github.com/danxhuber/isoclassify

# Download dependencies (see requirements.txt)

# Download MESA models into isoclassify directory
cd isoclassify
wget https://www.dropbox.com/s/dt1ts56gc9f7lzl/mesa.h5
wget https://www.dropbox.com/s/921jc0ojlz6c6ar/bcgrid.h5

# Set environment variables
export ISOCLASSIFY=yourpath/isoclassify # access mesa models via ${ISOCLASSIFY}/mesa.ebf 
export PYTHONPATH=yourpath/isoclassify:$PYTHONPATH
export PATH=yourpath/isoclassify/bin:$PATH # This adds isoclassify executable to your path
```

## Grid Modeling:

Derives posterior distributions for stellar parameters (teff, logg, feh, rad, mass, rho, lum, age, distance, av) through direct intregration of isochrones, given any set of observables (photometry, spectroscopy, parallax, asteroseismology) and priors. Can fit for extinction or use reddening map, includes (optionally) a correction to the Dnu scaling relation corrections by Sharma et al. (2016). <br />

See example/grid.ipynb for an application.

## Direct Method:

Uses bolometric corrections and extinction maps to derive stellar parameters using direct physical relations. Options are: <br />

1.  RA/DEC + Photometry + (Spectroscopy) + Parallax -> Teff, R, L, distance, Av. Uses Monte-Carlo method to implement exponentially decreasing volume density prior with a length scale of 1.35 kpc (Astraatmadja & Bailer-Jones 2016)

1. RA/DEC + Asteroseismology + Spectroscopy -> logg, rho, rad, mass, lum, distance, Av. Uses Dnu scaling relation corrections by Sharma et al. (2016) [not yet implemented]

Bolometric corrections are interpolated in (Teff, logg, FeH, Av) from the MIST grid (http://waps.cfa.harvard.edu/MIST/model_grids.html)


## Command line interface

isoclassify includes a command line interface (CLI) for convenient single star processing, as well as batch processing of many stars.

```bash
isoclassify run <mode> <star name> --csv <csv file> --outdir <output directory> --plot <plotmode> 
```

1. `<mode>` direct or grid
1. `<star name>` tells isoclassify which row of <csv file> to pull the input parameters
1. `<output directory>` will contain the isoclassify output files
1. `<csv file>` contains as columns parameters that are passed to isoclassify
1. `<plotmode>` tells isoclassify whether pop interactive plots `show`, save to a png `save-png`, or to not plot at all `none`.

We've included `examples/example.csv` with a few stars as an example

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
bin/isoclassify scrape-output 'output/*/output.csv' output.csv
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

## Testing the codebase

When making changes, if you haven't changed the core algorithm then run the following command

```
isoclassify batch direct examples/example.csv -o output/direct > isoclassify-test-direct.tot 
isoclassify batch grid examples/example.csv -o output/grid > isoclassify-test-grid.tot 
cat isoclassify-test-direct.tot isoclassify-test-grid.tot  > isoclassify-test.tot
parallel :::: isoclassify-test.tot 
```

and compare the output with a previous version.
