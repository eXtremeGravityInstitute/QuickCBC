import argparse
import glob
import matplotlib
matplotlib.use('Agg')
import corner
import numpy as np

data = np.loadtxt("spinpriors.dat", usecols=(0))

figure = corner.corner(data, labels=["$a$"], bins=40, title_kwargs={"fontsize": 14}, label_kwargs={"fontsize": 18}, show_titles=True)

figure.savefig("spinprior.png", dpi=300)

