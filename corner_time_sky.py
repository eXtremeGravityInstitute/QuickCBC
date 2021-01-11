import argparse
import glob
import matplotlib
matplotlib.use('Agg')
import corner
import numpy as np

data = np.loadtxt("allchain.dat", usecols=(6, 8, 9))

data[:,1] = data[:,1]*57.2957795130823209
data[:,2] = np.arccos(data[:,2])*57.2957795130823209

figure = corner.corner(data, labels=["$t_c$", "RA", "DEC"], bins=40, title_kwargs={"fontsize": 12}, label_kwargs={"fontsize": 14}, show_titles=True)

figure.savefig("skytime.png", dpi=300)
