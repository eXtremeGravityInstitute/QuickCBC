import argparse
import glob
import matplotlib
matplotlib.use('Agg')
import corner
import numpy as np

data = np.loadtxt("ratios.dat", usecols=(0, 1, 2))

figure = corner.corner(data, labels=[r'$\log_{10}(H/L)$', r'$\log_{10}(H/V)$', r'$\log_{10}(L/V)$'], bins=40, title_kwargs={"fontsize": 12}, label_kwargs={"fontsize": 14}, show_titles=True)

figure.savefig("corner_ratios.png", dpi=300)
