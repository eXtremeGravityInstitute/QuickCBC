#!/opt/local/bin/python2.7

import numpy as np
import healpy as hp
import pylab as pl
import matplotlib.pyplot as plt
from numpy import genfromtxt
from matplotlib import cm

grey_cmap = cm.Greys
grey_cmap.set_under("w") # sets background to white

orange_cmap = cm.Oranges
orange_cmap.set_under("w") # sets background to white

map = hp.read_map('sky.fits')
#Smoothes the map with a 1-degree FWHM Gaussian (fwhm given in radians).
map_smth = hp.sphtfunc.smoothing(map, fwhm = 0.017, iter = 1)


#hp.mollview(map_smth,title=" ", rot=[0,0,0], flip='astro', min=0, cbar=False, cmap=orange_cmap)
hp.mollview(map_smth,title=" ", flip='geo', min=0, cbar=False, cmap=orange_cmap)
hp.graticule(dpar=30,dmer=45)
#hp.graticule(dpar=15,dmer=30)

plt.savefig('sky.png')
