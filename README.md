# Isochrone Stellar Classification Codes

Python codes to perform stellar classifications given any set of input observables.

##Grid-Modeling:

Derives posterior distributions for stellar parameters through direct intregration of isochrones, given any set of observables (photometry, spectroscopy, parallax, asteroseismology) and priors. Can fit for extinction or use reddening map, includes (optionally) a correction to the Dnu scaling relation corrections by Sharma et al. (2016). <br />

Required isochrone grid calculated from MIST models (Choi et al. 2016): <br />
https://www.dropbox.com/s/yjgm8bwpw9cs0ew/mesa.ebf?dl=0 <br />

##Direct Method:

Ue MIST bolometric correction tables and extinction map to derive stellar parameters <br />
