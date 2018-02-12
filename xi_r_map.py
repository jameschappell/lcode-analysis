import h5py
import numpy as np
import matplotlib.pyplot as plt
import os

font = {'weight': 'normal', 'size': 15}
plt.rc('font', **font)

h5 = '/unix/pdpwa/jchappell/lcode/multi_beam_short/tb00010.h5'

def h5import(file):

    with h5py.File(file, 'r') as hf:

        print('Reading ' + file + '...')

        Q = hf.get('macroparticle_charge')[0]
        xi = np.array(hf.get('xi_mm'))
        x = np.array(hf.get('x_mm'))
        y = np.array(hf.get('y_mm'))
        px = np.array(hf.get('px_GeV'))
        py = np.array(hf.get('py_GeV'))
        pz = np.array(hf.get('pz_GeV'))

        print('Done.')
        print('Macroparticle charge = %.2e C' %Q)

    return Q, xi, x, y, px, py, pz


def plot_xi_r_map(xi_in, r_in, Q_in, filename, xi_min=-500.0, xi_max=0.0,
                  r_max=1.0, xi_bin=0.1, r_bin=0.01, vmax=3e-12, cmap='afmhot'):

    print('Estimating the 2D histogram...')
    H, xedges, yedges = np.histogram2d(xi_in, r_in,
                                       range=[[xi_min, xi_max], [0, r_max]],
                                       bins=[int((xi_max - xi_min) / xi_bin),
                                             int(r_max / r_bin)],
                                       weights = Q_in/(r_in+1e-10))

    # H needs to be rotated and flipped
    H = np.rot90(H)
    H = np.flipud(H)

    print('Plotting the 2D histogram...')
    fig = plt.figure(figsize=(10,3))

    plt.pcolormesh(xedges, yedges, H, cmap=cmap, vmin=0, vmax=vmax)
    plt.pcolormesh(xedges, -yedges, H, cmap=cmap, vmin=0, vmax=vmax)

    plt.xlabel(r'$\xi\ (\mathrm{mm})$')
    plt.ylabel(r'$x\ (\mathrm{mm})$')

    plt.ylim(-r_max, r_max)
    plt.xlim(xi_min, xi_max)

    plt.tight_layout()
    file_save_name = filename + '.png'
    plt.savefig(file_save_name)
    #plt.show()

    return H, xedges, yedges, xi_bin


# for saving beam images

numbers = np.arange(2000, 50001, 2000)
file_list = []
location = '/unix/pdpwa/jchappell/lcode/multi_beam/two_beam' \
              '/tb'

for i in numbers:
    number_string = '{0:05}'.format(i)
    string = location + number_string + '.h5'
    file_list.append(string)

for j in file_list:

    Q, xi, x, y, px, py, pz = h5import(j)
    r = np.sqrt(x * x + y * y)  # mm
    path_string = os.path.splitext(j)[0]
    H, xedges, yedges, xi_bin = plot_xi_r_map(xi, r, Q, path_string)

