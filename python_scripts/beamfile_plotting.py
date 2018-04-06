import matplotlib.pyplot as plt
import numpy as np
import os
import beamfile_output
import argparse

n = 7.0e14
c = 2.9979245800e8
w_p = 2 * np.pi * 8980 * np.sqrt(n)
E_0 = 0.511e-3 / (c/w_p)

#beamfileloc = '/unix/pdpwa/jchappell/lcode/multi_beam/plasma_state_1000/tb'

#files = ['10000']


def beamfile_plotting(location, lower, upper, step, xi_min, xi_max, xi_bin,
                      r_max, r_bin, batch):

    files = np.arange(lower, upper + 1, step)

    for i in files:

        number_string = '{0:05}'.format(i)
        beamfile = location + 'tb' + number_string + '.swp'

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

        #xi_min = -250.0
        #xi_max = -0.0
        #r_max = 1.0
        #r_bin = 0.01
        #xi_bin = 0.03

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

        if batch == 'False':

            plt.show()
            plt.savefig(figname)

        else:

            plt.savefig(figname)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""
            This script plots LCODE output beam files in hdf5 format as 
            .png files.""",
                                     formatter_class=
                                     argparse.RawTextHelpFormatter)

    parser.add_argument('--location', dest='location', default=None,
                        help='''
            This is the path to the directory of beam files that are to be 
            plotted.

            E.g. --file "<path to file>"''')

    parser.add_argument('--lower', dest='lower', default=None, help='''
            This is the value of the lowest beam file to loop over. For 
            example, 
            if your first beam file was 'tb00100.h5' then you would write:

            --lower 100
            ''')

    parser.add_argument('--upper', dest='upper', default=None, help='''
            This is the value of the highest beam file to loop over. For 
            example, 
            if your last beam file was 'tb20000.h5' then you would write:

            --upper 20000''')

    parser.add_argument('--step', dest='step', default=None, help='''
            This is the value of the step size between each beam file to loop 
            over. 
            For example, if each beam file was separated by 100 (i.e. 
            tb00100.h5, tb00200.h5, tb00300.h5 ...) then you would write:

            --step 100''')

    parser.add_argument('--xi_min', dest='xi_min', default=None, help='''
            This is the minimum value of xi to plot between.''')

    parser.add_argument('--xi_max', dest='xi_max', default=None, help='''
                This is the maximum value of xi to plot between.''')

    parser.add_argument('--xi_bin', dest='xi_bin', default=None, help='''
                This is the value of the desired bin width in xi for 
                plotting.''')

    parser.add_argument('--r_max', dest='r_max', default=None, help='''
                This is the maximum value of r to plot over.''')

    parser.add_argument('--r_bin', dest='r_bin', default=None, help='''
                This is the value of the desired bin width of r for 
                plotting.''')

    parser.add_argument('--batch', dest='batch', default=False, help='''
                This toggles batch mode for producing plots for a group of 
                images. Plot previews will not be shown. Can be either True 
                or False. Default is False.
                ''')

    arguments = parser.parse_args()

    # store input values

    location = arguments.location
    lower = int(arguments.lower)
    upper = int(arguments.upper)
    step = int(arguments.step)
    xi_min = float(arguments.xi_min)
    xi_max = float(arguments.xi_max)
    xi_bin = float(arguments.xi_bin)
    r_max = float(arguments.r_max)
    r_bin = float(arguments.r_bin)
    batch = arguments.batch

    # produce images

    beamfile_plotting(location, lower, upper, step, xi_min, xi_max, xi_bin,
                      r_max, r_bin, batch)