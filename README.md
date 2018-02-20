# Isochrone Stellar Classification Codes

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.573372.svg)](https://doi.org/10.5281/zenodo.573372)

Python codes to perform stellar classifications given any set of input observables. Details are described in [Huber et al. 2017](http://adsabs.harvard.edu/abs/2017ApJ...844..102H) - please cite this paper when using the code.

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
