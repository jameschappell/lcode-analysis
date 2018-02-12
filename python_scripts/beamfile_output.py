import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt


n = 7.0e14
c = 2.9979245800e8
w_p = 2 * np.pi * 8980 * np.sqrt(n)
E_0 = 0.511e-3 / (c/w_p)

print('c/w_p = {0:.2f} mm'.format(1000*c/w_p))
print('E_0 = {0:.2f} GV/m'.format(E_0))

L_plasma = 50000*c/w_p

beamfile = '../awake-baseline/beamfile.bin'

h5 = os.path.splitext(beamfile)[0] + '.h5'


def read_beamfile(bin, sample_fraction=1.0):

    dt = np.dtype([('xi', 'double'), ('r', 'double'), ('pz', 'double'),
                   ('pr', 'double'), ('M', 'double'), ('c2m', 'double'),
                   ('q', 'double'), ('id','double')])

    df = pd.DataFrame()
    fsize = os.path.getsize(bin)

    with open(bin, 'r') as fid:

        while True:

            data = np.fromfile(fid, dtype=dt, count=int(1e6))

            if len(data) <= 0:
                break

            print("\rReading beamfile ({:.0f} %)...".format(100 * fid.tell() /
                                                            fsize))
            df1 = pd.DataFrame(data.tolist(), columns=data.dtype.names)

            if sample_fraction < 1:
                df1 = df1.sample(frac=sample_fraction)

            df1[['id']] = df1[['id']].astype(long)
            df1 = df1.sort_values(by=['xi'], ascending=[False])

            df = df.append(df1, ignore_index=True)

    df = df.drop(df.tail(1).index)
    print('\nDone.')

    return df


def calculate_xy(df):

    # first calculate random angle:
    df['phi'] = 2 * np.pi * np.random.random(size=len(df))

    df['p_phi'] = df['M'] / (df['r'] + 1.0e-10)
    df['px'] = df['pr'] * np.cos(df['phi']) - df['p_phi'] * np.sin(df['phi'])
    df['py'] = df['pr'] * np.sin(df['phi']) + df['p_phi'] * np.cos(df['phi'])

    # converting to understandable units:
    mc = 0.511e-3   # GeV/c
    df['px_GeV'] = mc * df['px']
    df['py_GeV'] = mc * df['py']
    df['pz_GeV'] = mc * df['pz']
    df['xi_mm'] = df['xi'] * 1000 * (c/w_p)
    df['r_mm'] = df['r'] * 1000 * (c/w_p)
    df['x_mm'] = df['r_mm'] * np.cos(df['phi'])
    df['y_mm'] = df['r_mm'] * np.sin(df['phi'])

    df = df[['xi_mm', 'r_mm', 'x_mm', 'y_mm', 'px_GeV', 'py_GeV', 'pz_GeV',
             'q', 'c2m']]

    return df


def propagate(df, L=0):         # L -- propagation distance (in metres)

    df['x_mm'] = df['x_mm'] + L * 1e3 * df['px_GeV'] / df['pz_GeV']
    df['y_mm'] = df['y_mm'] + L * 1e3 * df['py_GeV'] / df['pz_GeV']
    df['r_mm'] = np.sqrt(df['x_mm'] * df['x_mm'] + df['y_mm'] * df['y_mm'])

    mc = 0.511e-3 / df['c2m']   # GeV/c
    p2 = df['px_GeV'] * df['px_GeV'] + df['py_GeV'] * df['py_GeV'] + \
         df['pz_GeV'] * df['pz_GeV']
    gamma = np.sqrt(1 + p2 / (mc*mc))
    beta_z = df['pz_GeV'] / (gamma*mc)

    df['xi_mm'] = df['xi_mm'] + (beta_z - 1) * L * 1e3

    return df


