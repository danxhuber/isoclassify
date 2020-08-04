import numpy as np
from astropy.io import ascii

out=ascii.read('output_direct.csv')

dust=out['id_starname'].astype(np.str)
dust[:]='allsky'
band=out['id_starname'].astype(np.str)
band[:]='vmag'
ascii.write([out['id_starname'],out['dir_teff'],out['dir_teff_err1'],out['dir_lum'],(np.abs(out['dir_lum_err1'])+np.abs(out['dir_lum_err2']))/2.,out['ra'],out['dec'],out['feh'],out['feh_err'],dust,band],'input_grid.csv',names=['id_starname','teff','teff_err','lum','lum_err','ra','dec','feh','feh_err','dust','band'],delimiter=',')

