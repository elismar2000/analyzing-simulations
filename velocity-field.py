import numpy as np
import matplotlib.pyplot as plt

from pygadgetreader import *

import sys


snapshot = sys.argv[1]

pixel_scale = np.deg2rad(0.2 / (60*60))
distance = 39e+3 #Kpc
physical_size = pixel_scale * distance

position = readsnap(snapshot, 'pos', 'gas')
velocity = readsnap(snapshot, 'vel', 'gas')

x = position[:, 0]
y = position[:, 1]

vz = velocity[:, 2]

# bins = (x.max() - x.min()) / physical_size



plt.hist2d(x, y, bins=100, weights=vz, cmap='Spectral')
plt.xlabel('x [Kpc]')
plt.ylabel('y [Kpc]')
cbar = plt.colorbar()
cbar.set_label('Vz [km/s]')
plt.show()
