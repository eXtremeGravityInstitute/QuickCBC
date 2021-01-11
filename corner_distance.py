import argparse
import glob
import matplotlib
matplotlib.use('Agg')
import corner
import numpy as np


data = np.loadtxt("allchain.dat", usecols=(11, 7))
#data = np.loadtxt("intrinsicchain.dat", usecols=(11, 7))

# convert to iota in degrees from cost(iota)
data[:,0] = np.arccos(data[:,0])*57.2957795130823209

figure = corner.corner(data, labels=["$\iota$","$D_L$"], bins=40, title_kwargs={"fontsize": 14}, label_kwargs={"fontsize": 16}, show_titles=True)

figure.savefig("distance.png", dpi=300)


