import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import unsio.input as uns


num_of_snaps = 44
path = '/home/elismar/Documentos/Fisica/IC/queorbita/orbits_9th_attempt/bigger_disk/'


def stars(snap_path):
    '''
    Read snapshot atributes for formed stars
    '''
    s = uns.CUNS_IN(snap_path,'all')
    s.nextFrame()

    _, ids_disk = s.getData('stars','id')

    if len(ids_disk) != 0:
        return True
    else:
        return False


disk_mass = []
star_mass = []
for snap in np.arange(0, num_of_snaps+1, 1):
    snapshot = 'snapshot_' + '%04d' % (snap,)
    snap_path = path + snapshot

    s = uns.CUNS_IN(snap_path,'all')
    s.nextFrame()

    _, mass_disk = s.getData('disk', 'mass')


    if stars(snap_path): _, mass_star = s.getData('stars', 'mass')
    else: mass_star = 0.0


    disk_mass.append(np.sum(mass_disk))
    star_mass.append(np.sum(mass_star))


time = np.arange(0, 45, 1) * 0.025 * 0.98

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

plt.style.use('seaborn-white')
cmap = plt.get_cmap('gist_stern')
colors = [cmap(i) for i in np.linspace(0, 1, 7)]

fig = plt.figure()
axs = fig.add_subplot()

axs.plot(time, disk_mass, color=colors[1], label='Mass in the disk')
axs.plot(time, star_mass, color=colors[2], label='Mass in recently formed stars')
axs.set_xlabel('Time [Gyr]')
axs.set_ylabel(r'Stellar Mass [$10^{10}\ M_{\odot}$]')
axs.axvline(1.0, color='red', label=r'1.0Gyr $\approx$ perigalacticon')
axs.set_title('NGC2992')
plt.legend()
plt.show()
