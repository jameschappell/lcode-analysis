import numpy as np
import pandas as pd
import os
import h5py
import holoviews as hv
from IPython.display import display_png

import matplotlib.pyplot as plt
import matplotlib


bin = '/Users/jamiechappell/Documents/PhD/Simulations/lcode/PBC/dephasing' \
      '/parameter_scans/charge/1000pC/beamfile.bin'
t = 1000000
sample_fraction = 1.0

dt = np.dtype(
    [('xi', 'double'), ('r', 'double'), ('pz', 'double'), ('pr', 'double'),
     ('M', 'double'), ('c2m', 'double'), ('q', 'double'), ('id', 'double')])

df = pd.DataFrame()
fsize = os.path.getsize(bin)
with open(bin, 'r') as fid:
    for i in range(1, 15):
        data = np.fromfile(fid, dtype=dt, count=-1)
        if len(data) > 0:
            print("\rReading beamfile.bin ({:.0f} %)...".format(
                100 * fid.tell() / fsize))
            df1 = pd.DataFrame(data.tolist(), columns=data.dtype.names)
            if sample_fraction < 1:
                df1 = df1.sample(frac=sample_fraction)
            df1[['id']] = df1[['id']].astype(long)
            df1 = df1.sort_values(by=['xi'], ascending=[False])

            df = df.append(df1, ignore_index=True)

print('\nDone.')

df = df.drop(df.tail(1).index)

df = df[df.c2m > 0.5]

n = 2.8e15
c = 2.9979245800e8
w_p =2*np.pi*8980*np.sqrt(n)
E_0 = 0.511e-3/(c/w_p)
print('c/w_p = {0:.2f} mm'.format(1000*c/w_p))
print('E_0 = {0:.2f} GV/m'.format(E_0))

xi = (df['xi']*1e3*c/w_p).values # mm
pz = (df['pz']*0.511e-3).values  # GeV

dim_xi = hv.Dimension('xi', unit='mm', range=(-1.0,0.0))
dim_pz = hv.Dimension('pz', unit='GeV', range=(0,100.0))
#xi_r  = hv.Points(df[['xi_mm', 'r_mm' ]])

print df

matplotlib.get_backend()


plt.scatter(xi, pz)
plt.show()

