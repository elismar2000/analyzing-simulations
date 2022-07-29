import numpy as np
import matplotlib.pyplot as plt
import unsio.input as uns

path = '/home/elismar/Documentos/Fisica/IC/queorbita/orbits_13th_attempt/orb17-fric/'
snap = 42


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


def _distances(position, point):
    '''
    Evaluates distances of particles from a given point
    '''

    distances = np.sqrt(np.sum(np.square(position[i] - point[i]) for i in range(3)))
    return distances


def calculate_Re(path, snap):


    snapshot = 'snapshot_' + '%04d' % (snap,)
    snap_path = path + snapshot
    print("Calculating for snapshot ", snapshot)

    data = read_basic(snap_path)
    pos = data['position']
    ids = data['ids']
    mass = data['mass']

    mask2992 = (ids >= 350000) & (ids <= 362499)
    mask2993 = (ids >= 362500) & (ids <= 374999)

    disk_mass_2992 = np.sum(mass[mask2992])
    disk_mass_2993 = np.sum(mass[mask2993])

    x = pos[0::3]
    y = pos[1::3]
    z = pos[2::3]

    pos = np.vstack((x, y, z)).T

    coords_min2992, coords_min2993 = _coords_min(snap_path, case='interaction')

    dists2992 = np.array([_distances(pos[mask2992][i], coords_min2992) for i in range(len(pos[mask2992]))])
    dists2993 = np.array([_distances(pos[mask2993][i], coords_min2993) for i in range(len(pos[mask2993]))])

    Re2992 = []
    Re2993 = []
    #Raio efetivo = metade da massa
    for r in np.linspace(3.0, 8.0, 300):
        mass_inside = np.sum(mass[mask2992][dists2992 < r])

        if (mass_inside > (0.264 * disk_mass_2992)):
            Re2992.append(r)
            break

    import pdb; pdb.set_trace()
    for r in np.linspace(2.0, 8.0, 300):
        mass_inside = np.sum(mass[mask2993][dists2993 < r])

        if (mass_inside > (0.264 * disk_mass_2993)):
            Re2993.append(r)
            break

    return Re2992, Re2993


snapshot = 'snapshot_' + '%04d' % (snap,)
snap_path = path + snapshot
data = read_basic(snap_path)

coords_min2992, coords_min2993 = _coords_min(snap_path, 'interaction')

pos = data['position']
ids = data['ids']

mask2992 = (ids >= 350000) & (ids <= 362499)
mask2993 = (ids >= 362500) & (ids <= 374999)

x = pos[0::3]
y = pos[1::3]
z = pos[2::3]

fig, axs = plt.subplots(2, 1)
axs[0].plot(x[mask2992], y[mask2992], ',')
axs[0].plot(coords_min2992[0], coords_min2992[1], 'rX')
axs[0].set_xlabel('x')
axs[0].set_ylabel('y')

axs[1].plot(x[mask2992], z[mask2992], ',')
axs[1].set_xlabel('x')
axs[1].set_ylabel('z')
plt.show()
