import h5py
import beamfile_output
import numpy as np
import os
import argparse


n = 7.0e14
c = 2.9979245800e8
w_p = 2 * np.pi * 8980 * np.sqrt(n)
E_0 = 0.511e-3 / (c/w_p)

#numbers = np.arange(34000, 50001, 2000)

#location = '/unix/pdpwa/jchappell/lcode/multi_beam/two_beam/tb'


def save_to_hdf5(location, lower, upper, step):

    file_list = []

    numbers = np.arange(lower, upper, step)

    for i in numbers:
        number_string = '{0:05}'.format(i)
        string = location + 'tb' + number_string + '.swp'
        file_list.append(string)


    for beamfileloc in file_list:

        h5 = os.path.splitext(beamfileloc)[0] + '.h5'

        df0 = beamfile_output.read_beamfile(beamfileloc)
        df = beamfile_output.calculate_xy(df0)

        dxi = 0.01 * (c/w_p)        # m (simulation step code used in lcode.cfg)
        re = 2.818e-15              # m
        q0 = 1.6e-19                # C
        Q = ((dxi/re) * q0 / 2) * df['q'].head(1).values
        # macroparticle charge Q

        print("Saving beamfile to " + h5 + "...")
        os.path.exists(h5) and os.remove(h5)

        with h5py.File(h5, 'w') as hf:
            hf.create_dataset("macroparticle_charge", data=Q, dtype='float32')
            # in C

            for col in ['x_mm', 'y_mm', 'xi_mm', 'px_GeV', 'py_GeV', 'pz_GeV']:

                print(col+' ')
                ds = hf.create_dataset(col, data=df[col].values,
                                       dtype='float32')

    print("\nDone.")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""
        This script converts LCODE output beam files from swp files to hdf5 
        files for quicker processing.""",
                                     formatter_class=
                                     argparse.RawTextHelpFormatter)

    parser.add_argument('--location', dest='location', default=None,
                        help='''
        This is the path to the directory of beam files that are to be 
        converted.

        E.g. --file "<path to file>"''')

    parser.add_argument('--lower', dest='lower', default=None, help='''
        This is the value of the lowest beam file to loop over. For example, 
        if your first beam file was 'tb00100.swp' then you would write:
        
        --lower 100
        ''')

    parser.add_argument('--upper', dest='upper', default=None, help='''
        This is the value of the highest beam file to loop over. For example, 
        if your last beam file was 'tb20000.swp' then you would write:
        
        --upper 20000''')

    parser.add_argument('--step', dest='step', default=None, help='''
        This is the value of the step size between each beam file to loop over. 
        For example, if each beam file was separated by 100 (i.e. 
        tb00100.swp, tb00200.swp, tb00300.swp...) then you would write:
        
        --step 100''')

    arguments = parser.parse_args()

    # store values given in input

    location = arguments.location
    lower = int(arguments.lower)
    upper = int(arguments.upper)
    step = int(arguments.step)

    # convert beam files

    save_to_hdf5(location, lower, upper, step)


