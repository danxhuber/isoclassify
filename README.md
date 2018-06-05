# Isochrone Stellar Classification Codes

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.573372.svg)](https://doi.org/10.5281/zenodo.573372)

Python codes to perform stellar classifications given any set of input observables. Details are described in [Huber et al. 2017](http://adsabs.harvard.edu/abs/2017ApJ...844..102H) - please cite this paper when using the code.

## Installation

```bash
# Clone repo
git clone https://github.com/danxhuber/isoclassify

# Download dependencies
pip install ebfpy 
pip install ephem
pip install dustmaps

# Download MESA models into isoclassify directory
cd isoclassify
wget https://www.dropbox.com/s/yjgm8bwpw9cs0ew/mesa.ebf?dl=0 

# Download asfgrid input files for asteroseismic Delta_nu corrections
cd isoclassify/direct
wget https://www.dropbox.com/s/26enigp5q03bz5n/grid_interp1.ebf?dl=0
wget https://www.dropbox.com/s/xrwwhgdgrhwgms7/grid_interp2.ebf?dl=0

# Set environment variables
export ISOCLASSIFY=${WKDIR}/code/isoclassify # access mesa models via ${ISOCLASSIFY}/mesa.ebf 
export PYTHONPATH=${PYTHONPATH}:${ISOCLASSIFY}
export PATH=${WKDIR}/code/isoclassify/bin:$PATH # This adds isoclassify executable to your path
```

## Grid Modeling:

Derives posterior distributions for stellar parameters (teff, logg, feh, rad, mass, rho, lum, age, distance, av) through direct intregration of isochrones, given any set of observables (photometry, spectroscopy, parallax, asteroseismology) and priors. Can fit for extinction or use reddening map, includes (optionally) a correction to the Dnu scaling relation corrections by Sharma et al. (2016). <br />

Required isochrone grid interpolated from MIST models (Choi et al. 2016): <br />
https://www.dropbox.com/s/yjgm8bwpw9cs0ew/mesa.ebf?dl=0 <br />

See grid/example.ipynb for an application.

## Direct Method:

Uses bolometric corrections and extinction maps to derive stellar parameters using direct physical relations. Options are: <br />
(1) RA/DEC + Photometry + (Spectroscopy) + Parallax -> Teff, R, L, distance, Av. Uses Monte-Carlo method to implement exponentially decreasing volume density prior with a length scale of 1.35 kpc (Astraatmadja & Bailer-Jones 2016) <br />
(2) RA/DEC + Asteroseismology + Spectroscopy -> logg, rho, rad, mass, lum, distance, Av. Uses Dnu scaling relation corrections by Sharma et al. (2016). <br />

Bolometric corrections are interpolated in (Teff, logg, FeH, Av) from the MIST grid, Conroy et al. in prep (http://waps.cfa.harvard.edu/MIST/model_grids.html)


## Command line interface

isoclassify includes a command line interface (CLI) for convenient single star processing, as well as batch processing of many stars.

```bash
isoclassify run <mode> <star name> --csv <csv file> --outdir <output directory>  
```

1. `<mode>` direct or grid
1. `<star name>` tells isoclassify which row of <csv file> to pull the input parameters
1. `<output directory>` will contain the isoclassify output files
1. `<csv file>` contains as columns parameters that are passed to isoclassify

We've included `examples/example.csv` with a few stars as an example

```bash
mkdir -p output/sol # make sure output directory exists
isoclassify run direct sol --csv examples/example.csv --outdir output/sol
```

The CLI also makes parallel processing easy. First generate a list of batch scripts (isoclassify.tot), make isoclassify.tot an executable, then run using GNU parallel

```bash
isoclassify batch direct examples/example.csv -o output > isoclassify.tot
chmod +x isoclassify.tot
parallel ::: isoclassify.tot
```

