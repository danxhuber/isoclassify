import numpy as np
from astropy.io import ascii

out=ascii.read('output_grid.csv')
inp=ascii.read('input_direct.csv')

inp['logg']=out['iso_logg']
inp['logg_err']=out['iso_logg_err1']

ascii.write(inp,'input_direct.csv',delimiter=',')