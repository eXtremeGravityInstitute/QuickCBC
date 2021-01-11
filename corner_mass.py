import argparse
import glob
import matplotlib
matplotlib.use('Agg')
import corner
import numpy as np

data = np.loadtxt("allchain.dat", usecols=(15, 16, 14, 13, 12, 7))

figure = corner.corner(data, labels=[r'$m_1$', r'$m_2$', r'$M$', r'${\cal M}$', r'$z$', r'$D_L$'], bins=40, title_kwargs={"fontsize": 14}, label_kwargs={"fontsize": 18}, show_titles=True)

figure.savefig("masses.png", dpi=300)

