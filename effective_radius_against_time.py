import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import unsio.input as uns
from tables import *


def read_basic(snapshot):
    """
    Reads snapshots and puts them in nice dictionaries.
    Parameters
    ----------
    snapshot : str
        Name of the file containing the snapshot.

    Returns
    -------
    data : dict
        Dictionary containing all the relevant information
        from the snapshot.
    """

    s = uns.CUNS_IN(snapshot,'all')
    s.nextFrame()

    ok, pos  = s.getData('disk','pos')
    ok, ids  = s.getData('disk','id')
    ok, mass  = s.getData('disk','mass')
    ok, pot = s.getData('disk', 'pot')

    data = {'position': pos, 'ids': ids, 'mass': mass, 'pot': pot,}

    return data


def _distances(position, point):
    '''
    Evaluates distances of particles from a given point
    '''

    distances = np.sqrt(np.sum(np.square(position[i] - point[i]) for i in range(3)))
    return distances


def _coords_min(snapshot, case):
    '''
    Calculate coordinates of the minimum of potential
    of a galaxy
    '''

    s = uns.CUNS_IN(snapshot,'all')
    s.nextFrame()

    data = read_basic(snapshot)
    ids = data['ids']
    pos = data['position']
    pot = data['pot']


    if case == 'isolated':
        pot1 = pot.min()
        index_min = np.where(pot == pot1)[0][0]
        coords_min = pos[index_min*3:index_min*3+3]

        return coords_min

    if case == 'interaction':
        mask2992 = (ids >= 350000) & (ids <= 362499)
        mask2993 = (ids >= 362500) & (ids <= 374999)

        pot_min2992 = pot[mask2992].min()
        index_min2992 = np.where(pot == pot_min2992)[0][0]
        coords_min2992 = pos[index_min2992*3:index_min2992*3+3]

        pot_min2993 = pot[mask2993].min()
        index_min2993 = np.where(pot == pot_min2993)[0][0]
        coords_min2993 = pos[index_min2993*3:index_min2993*3+3]

        return coords_min2992, coords_min2993


def calculate_Re(num_of_snaps, path, orbit):

    Re2992 = np.zeros(45)
    Re2993 = np.zeros(45)

    for snap in range(0, num_of_snaps+1, 1):
        snapshot = 'snapshot_' + '%04d' % (snap,)
        snap_path = path + snapshot
        print("Calculating for snapshot ", snapshot)


        #Read basic properties
        data = read_basic(snap_path)
        pos = data['position']
        ids = data['ids']
        mass = data['mass']


        #Define masks for the disks of each galaxy
        mask2992 = (ids >= 350000) & (ids <= 362499)
        mask2993 = (ids >= 362500) & (ids <= 374999)


        #Calculate total mass of the disks
        disk_mass_2992 = np.sum(mass[mask2992])
        disk_mass_2993 = np.sum(mass[mask2993])


        #Write position vector in the right format
        x = pos[0::3]
        y = pos[1::3]
        z = pos[2::3]
        pos = np.vstack((x, y, z)).T


        #Calculate minima of potential of each galaxy disk
        coords_min2992, coords_min2993 = _coords_min(snap_path, case='interaction')

        if (snap == 39):
            if orbit == 3: coords_min2992 = np.array([8.125, 5.577, 3.119])
            if orbit == 13: coords_min2992 = np.array([8.735, 5.136, 4.022])
            if orbit == 17: coords_min2992 = np.array([7.279, 6.084, 3.967])


        #Calculate distances of each particle to the disk center
        dists2992 = np.array([_distances(pos[mask2992][i], coords_min2992) for i in range(len(pos[mask2992]))])
        dists2993 = np.array([_distances(pos[mask2993][i], coords_min2993) for i in range(len(pos[mask2993]))])


        #Effective radius = radius containing half the total mass
        for r in np.linspace(3.0, 8.0, 300):
            mass_inside = np.sum(mass[mask2992][dists2992 < r])

            if (mass_inside > (0.264 * disk_mass_2992)):
                Re2992[snap] = r
                break

        for r in np.linspace(2.0, 8.0, 300):
            mass_inside = np.sum(mass[mask2993][dists2993 < r])

            if (mass_inside > (0.264 * disk_mass_2993)):
                Re2993[snap] = r
                break

    return Re2992, Re2993


def calculate_Re_isolated(path, num_of_snaps):
    '''
    Calculate Re for isolated galaxies
    '''

    Re = []
    for snap in range(0, num_of_snaps+1, 1):
        snapshot = 'snapshot_' + '%04d' % (snap,)
        snap_path = path + snapshot
        print("Calculating for snapshot ", snapshot)

        #Read basic properties
        data = read_basic(snap_path)
        pos = data['position']
        mass = data['mass']


        #Write position vector in the right format
        x = pos[0::3]
        y = pos[1::3]
        z = pos[2::3]
        pos = np.vstack((x, y, z)).T


        #Calculate total mass of the disk
        disk_mass = np.sum(mass)
        # print('total disk mass = ', disk_mass)

        #Calculate position of the minimum of potential
        coords_min = _coords_min(snap_path, case='isolated')


        #Calculate distances to the minimum of potential
        dists = np.array([_distances(pos[i], coords_min) for i in range(len(pos))])


        #Calculate effetive radius Re
        for r in np.linspace(2.0, 6.0, 250):
           mass_inside = np.sum(mass[dists < r])

           if (mass_inside > (0.264 * disk_mass)):
               Re.append(r)
               break

    return np.asarray(Re)


def plot():

    ngc2992 = '/home/elismar/Documentos/Fisica/IC/Gadget3/testing_isolated_galaxies/ngc2992_MakeNewDisk/'
    ngc2993 = '/home/elismar/Documentos/Fisica/IC/Gadget3/testing_isolated_galaxies/ngc2993_MakeNewDisk/'
    num_of_snaps = 44

    Re2992_isolated = calculate_Re_isolated(ngc2992, num_of_snaps)
    Re2993_isolated = calculate_Re_isolated(ngc2993, num_of_snaps)

    time = np.arange(0, 45, 1) * 0.025 * 0.98

    from matplotlib import rc
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)


    from matplotlib.pyplot import cm
    colors = cm.gnuplot(np.linspace(0, 1, 3))
    # colors2 = cm.autumn(np.linspace(0, 1, 3))

    fig, axs = plt.subplots(2, 1)

    for orbit, c in zip([3, 13, 17], colors):
        path = '/home/elismar/Documentos/Fisica/IC/queorbita/orbits_13th_attempt/orb' + str(orbit) + '-fric/'
        print('This is orbit ', orbit)
        Re2992, Re2993 = calculate_Re(num_of_snaps, path, orbit)

        axs[0].plot(time, Re2992/Re2992_isolated, color=c)
        axs[1].plot(time, Re2993/Re2993_isolated, color=c,label='Orbit ' + str(orbit))


    axs[0].set_xlabel('Time [Gyr]')
    axs[0].set_ylabel(r'$R_{e,int} / R_{e,iso}$')
    axs[0].axvline(1.0)
    axs[0].axhline(1.0, linestyle='--')
    axs[0].set_title('Re2992-interaction/NGC2992-isolated')

    axs[1].set_xlabel('Time [Gyr]')
    axs[1].set_ylabel(r'$R_{e,int} / R_{e,iso}$')
    axs[1].axvline(1.0)
    axs[1].axhline(1.0, linestyle='--')
    axs[1].set_title('Re2993-interaction/NGC2993-isolated')
    axs[1].legend()

    plt.subplots_adjust(wspace=0.0, hspace=0.5)
    plt.show()



if __name__ == '__main__':

    plot()
