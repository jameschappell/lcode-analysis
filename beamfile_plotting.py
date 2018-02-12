import matplotlib.pyplot as plt
import numpy as np
import os
import beamfile_output

n = 7.0e14
c = 2.9979245800e8
w_p = 2 * np.pi * 8980 * np.sqrt(n)
E_0 = 0.511e-3 / (c/w_p)

beamfileloc = '/unix/pdpwa/jchappell/lcode/multi_beam/plasma_state_1000/tb'

files = ['10000']

for i in files:

    beamfile = beamfileloc + i + '.swp'

    h5 = os.path.splitext(beamfile)[0] + '.h5'
    name = os.path.splitext(beamfile)[0]
    print name
    df0 = beamfile_output.read_beamfile(beamfile)
    df = beamfile_output.calculate_xy(df0)

    font = {'weight': 'normal', 'size': 15}
    plt.rc('font', **font)

    uppath = lambda _path, n: os.sep.join(_path.split(os.sep)[n:])

    xi = df.xi_mm.values
    r = df.r_mm.values

    dxi = 0.01 * (c/w_p)        # m (simulation step code used in lcode.cfg)
    re = 2.818e-15              # m
    q0 = 1.6e-19                # C
    Q = ((dxi/re) * q0 / 2) * df['q'].head(1).values    # macroparticle charge Q

    xi_min = -250.0
    xi_max = -0.0

    r_max = 1.0


    r_bin = 0.01
    xi_bin = 0.03

    print('Estimating the 2D histogram...')
    H, xedges, yedges = np.histogram2d(xi, r, range=[[xi_min, xi_max],
                                                     [0, r_max]],
                                       bins=[int((xi_max - xi_min) / xi_bin),
                                             int(r_max / r_bin)],
                                       weights=Q/(r + 0.00001))

    # H needs to be rotated and flipped
    H = np.rot90(H)
    H = np.flipud(H)

    H_max = H.max()

    H = H / H_max

    print('Plotting the 2D histogram...')
    fig = plt.figure(figsize=(16, 5))

    plt.pcolormesh(xedges, yedges, H, cmap='hot', vmin=0, vmax=0.1)
    plt.pcolormesh(xedges, -yedges, H, cmap='hot', vmin=0, vmax=0.1)
    plt.xlabel(r'$\xi\ (\mathrm{mm})$')
    plt.ylabel(r'$x\ (\mathrm{mm})$')

    title = uppath(h5, 8)
    plt.title(title)
    plt.ylim(-r_max, r_max)
    plt.xlim(xi_min, xi_max)

    figname = name + '.png'
    plt.show()
    #plt.savefig(figname)


