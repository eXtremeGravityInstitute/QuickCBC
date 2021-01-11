import argparse
import glob
import matplotlib
matplotlib.use('Agg')
import corner
import numpy as np

data = np.loadtxt("allchain.dat", usecols=(8, 9))

figure = corner.corner(data, labels=["RA", "sin(DEC)"], bins=40, title_kwargs={"fontsize": 12}, label_kwargs={"fontsize": 14}, show_titles=True)

figure.savefig("cornersky.png", dpi=300)
