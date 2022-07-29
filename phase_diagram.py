import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import unsio.input as uns_in


snapshot = '/home/elismar/Documentos/Fisica/IC/Gadget3/testing_isolated_galaxies/dice_standard_nowarp_nospiral/snapshot_0000'

#=====================================
#Define functions to calculate r and p
#=====================================
def r(x, y, z):
    return np.sqrt(x**2 + y**2 + z**2)

def p(vx, vy, vz, m):
    # return m * np.sqrt(vx**2 + vy**2 + vz**2)
    return np.sqrt(vx**2 + vy**2 + vz**2)

#=====================================
#Read positions, velocities and masses of the particles
#=====================================
s = uns_in.CUNS_IN(snapshot,'all')
s.nextFrame()

def read_properties(component: str):
    _, pos = s.getData(component, 'pos')
    _, vel = s.getData(component, 'vel')
    _, mass = s.getData(component, 'mass')

    x = pos[0::3]
    y = pos[1::3]
    z = pos[2::3]

    vx = vel[0::3]
    vy = vel[1::3]
    vz = vel[2::3]

    return x, y, z, vx, vy, vz, mass


x_halo, y_halo, z_halo, vx_halo, vy_halo, vz_halo, mass_halo = read_properties('halo')
x_disk, y_disk, z_disk, vx_disk, vy_disk, vz_disk, mass_disk = read_properties('disk')
x_gas, y_gas, z_gas, vx_gas, vy_gas, vz_gas, mass_gas = read_properties('gas')
x_bulge, y_bulge, z_bulge, vx_bulge, vy_bulge, vz_bulge, mass_bulge = read_properties('bulge')


r_halo = r(x_halo, y_halo, z_halo)
p_halo = p(vx_halo, vy_halo, vz_halo, mass_halo)

r_disk = r(x_disk, y_disk, z_disk)
p_disk = p(vx_disk, vy_disk, vz_disk, mass_disk)

r_gas = r(x_gas, y_gas, z_gas)
p_gas = p(vx_gas, vy_gas, vz_gas, mass_gas)

r_bulge = r(x_bulge, y_bulge, z_bulge)
p_bulge = p(vx_bulge, vy_bulge, vz_bulge, mass_bulge)

#=====================================
#Plot
#=====================================

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

plt.style.use('seaborn-white')

cmap = plt.get_cmap('Set1')
colors = [cmap(i) for i in np.linspace(0, 1, 8)]

plt.scatter(r_halo, p_halo, s=2.0, color=colors[0], label='Halo', alpha=0.5)
plt.scatter(r_disk, p_disk, s=2.0, color=colors[2], label='Disk', alpha=0.5)
plt.scatter(r_gas, p_gas, s=2.0, color=colors[4], label='Gas', alpha=0.5)
plt.scatter(r_bulge, p_bulge, s=2.0, color=colors[6], label='Bulge', alpha=0.5)

plt.title('Dice standard (no warps, no spirals)')
plt.xlabel('Position [Kpc]')
# plt.ylabel(r'Momentum [$10^{10}\ M_{\odot}\ km\ s^{-1}$]')
plt.ylabel(r'Velocity [$km\ s^{-1}$]')

plt.xlim(0, 50.0)

plt.legend()
plt.show()
