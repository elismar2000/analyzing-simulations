import numpy as np
import matplotlib.pyplot as plt
import pygadgetreader as gadred
from astropy.io import fits
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel


def density_map(posxy, gal_prop, x_c, y_c, method='mean', scpix=1, ima_size=[150,150]):

    gal_den=np.zeros((ima_size[1],ima_size[0]))*np.nan

    x_c = (ima_size[0]*0.5)*scpix-x_c
    y_c = (ima_size[1]*0.5)*scpix-y_c

    for ind, e in np.ndenumerate(gal_den):
        cond = ((posxy[:,1] > ind[0]*scpix-y_c)  & (posxy[:,1] < ind[0]*scpix-y_c+scpix) &
                (posxy[:,0] > ind[1]*scpix-x_c)  & (posxy[:,0] < ind[1]*scpix-x_c +scpix) )
        if method == 'mean':
            gal_den[ind] = np.nanmean(gal_prop[cond])
        else:
            gal_den[ind] = np.nansum(gal_prop[cond])

    if method == 'sum':
        gal_den[gal_den==0.0]=np.nan

    return gal_den

snapshot = '../orbits_16th_attempt/orb147-fric-e0.9/snapshot_021'
gadred.readheader(snapshot,'header')

gas_vel = gadred.readsnap(snapshot,'vel','gas')
gas_pos = gadred.readsnap(snapshot,'pos','gas')
gas_mass = gadred.readsnap(snapshot,'mass','gas')
disk_pos = gadred.readsnap(snapshot,'pos','disk')
disk_mass = gadred.readsnap(snapshot,'mass','disk')

#===============================
#Plotting stars
#===============================
# hb = plt.hexbin(disk_pos[:,0],disk_pos[:,1], gridsize=50, bins='log', cmap='inferno')
# cb = plt.colorbar(hb)
# cb.set_label('log10(N)')
#
# #===============================
# #Plotting gas
# #===============================
# hb = plt.hexbin(gas_pos[:,0],gas_pos[:,1], gridsize=500, bins='log', cmap='inferno')
# cb = plt.colorbar(hb)
# cb.set_label('log10(N)')
# limx = np.percentile(disk_pos[:,0], [0,100])
# limy = np.percentile(disk_pos[:,1], [0,100])
# plt.xlim(limx)
# plt.ylim(limy)
#
# #===============================
# #Determining the center of the image
# #===============================
# xcent = -10
# ycent = -20
#
# plt.vlines(xcent,limy[0], limy[1])
# plt.hlines(ycent,limx[0], limx[1])
#
# hb = plt.hexbin(disk_pos[:,0],disk_pos[:,1], gridsize=50, bins='log', cmap='inferno')
# cb = plt.colorbar(hb)
# cb.set_label('log10(N)')
# plt.xlim(limx)
# plt.ylim(limy)
#
# #===============================
# #Making Star and gas images
# #===============================
dimx = 300
dimy = 300
print('here')
map_rad_vel = density_map(gas_pos, gas_vel[:,2] , -10, -20, ima_size=[dimx, dimy], scpix=0.5)
map_star=density_map(disk_pos, disk_mass, -10, -20, 'sum', ima_size=[dimx, dimy], scpix=0.5)
map_gas=density_map(gas_pos, gas_mass, -10, -20, 'sum', ima_size=[75, 75], scpix=2)

#===============================
#Making Star and gas images
#===============================
plt.figure(figsize=(7,7))
cb=plt.imshow(map_rad_vel, origin='lower', cmap='Spectral', vmin=-600, vmax=300, interpolation='none',
             extent=[-dimx,dimx,-dimy,dimy])

cb = plt.colorbar(cb)

plt.imshow(map_star, origin='lower', vmin=np.nanpercentile(map_star,10), alpha=0.5,
                                     vmax=np.nanpercentile(map_star,99.8), interpolation='none',
           extent=[-dimx,dimx,-dimy,dimy])



plt.minorticks_on()
plt.xlabel("Delta x (kpc)")
plt.ylabel("Delta y (kpc)")

plt.show()
