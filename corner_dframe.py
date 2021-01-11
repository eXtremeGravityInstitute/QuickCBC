import argparse
import glob
import matplotlib
matplotlib.use('Agg')
import corner
import numpy as np

data = np.loadtxt("allchain.dat", usecols=(2, 3, 19, 4))

figure = corner.corner(data, labels=[r'${\cal M}_{d}$', r'$M_{d}$', r'$t_{0}$', r'$\chi_{eff}$'], bins=40, title_kwargs={"fontsize": 12}, label_kwargs={"fontsize": 14}, show_titles=True)

figure.savefig("dframe.png", dpi=300)
