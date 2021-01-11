import argparse
import glob
import matplotlib
matplotlib.use('Agg')
import corner
import numpy as np

data = np.loadtxt("allchain.dat", usecols=(15, 16, 14, 4, 12, 7))

# get mass ratio to q = m2/m1
data[:,2] = data[:,1]/data[:,0]

figure = corner.corner(data, labels=["$m_1$", "$m_2$", "$q$", "$\chi_{eff}$", "$z$", "$D_L$"], bins=40, title_kwargs={"fontsize": 14}, label_kwargs={"fontsize": 18}, show_titles=True)

figure.savefig("mass_chi.png", dpi=300)

