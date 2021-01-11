import argparse
import glob
import matplotlib
matplotlib.use('Agg')
import corner
import numpy as np

data = np.loadtxt("allchain.dat", usecols=(8, 9, 19, 20))

data[:,0] = data[:,0]*57.2957795130823209
data[:,1] = np.arccos(data[:,1])*57.2957795130823209
data[:,2] = (data[:,2]-data[:,3])*1000.0

dx = np.delete(data, 3, axis=1)

figure = corner.corner(dx, labels=["RA", "DEC", "dt"], bins=40, title_kwargs={"fontsize": 12}, label_kwargs={"fontsize": 14}, show_titles=True)

figure.savefig("skytime.png", dpi=300)
