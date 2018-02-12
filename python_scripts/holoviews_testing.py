import numpy as np
import pandas as pd
import holoviews as hv
import os
from holoviews.operation.datashader import aggregate, shade, datashade, \
    dynspread
from holoviews.plotting.mpl import MPLRenderer

hv.notebook_extension('bokeh')

bin = '../awake-baseline/beamfile.bin'

sample_fraction = 1.0

dt = np.dtype([('xi', 'double'), ('r', 'double'), ('pz', 'double'),
               ('pr', 'double'), ('M', 'double'), ('c2m', 'double'),
               ('q', 'double'), ('id', 'double')])

df = pd.DataFrame()
fsize = os.path.getsize(bin)

with open(bin, 'r') as fid:

    for i in range(1,15):

        data = np.fromfile(fid, dtype=dt, count=int(1e6))

        if len(data) > 0:

            print("\rReading beamfile.bin ({:.0f} %)...".format(100*fid.tell() /
                                                                fsize))
            df1 = pd.DataFrame(data.tolist(), columns=data.dtype.names)

            if sample_fraction < 1:
                df1.sample(frac=sample_fraction)

            df1[['id']] = df1[['id']].astype(long)
            df1 = df1.sort_values(by=['xi'], ascending=[False])

            df = df.append(df1, ignore_index=True)


print('\nDone.')

df['pr_MeV'] = df['pr']*0.511
df['pz_MeV'] = df['pz']*0.511

xi_r = hv.Points(df[['xi', 'r']])
xi_pz = hv.Points(df[['xi', 'pz_MeV']])
xi_pr = hv.Points(df[['xi', 'pr_MeV']])

plot_1 = (dynspread(datashade(xi_r)) + dynspread(datashade(xi_pz)) +
 dynspread(datashade(xi_pr))).cols(1)
renderer = hv.renderer('matplotlib').instance(fig='svg', holomap='gif')
renderer.save(plot_1, 'example_I', style=dict(Image={'cmap':'jet'}))

# pr_r = hv.Points(df[['pr_MeV', 'r']])
# pr_pz = hv.Points(df[['pr_MeV', 'pz_MeV']])
#
# (dynspread(datashade(xi_r)) + dynspread(datashade(pr_r))).cols(2)
# (dynspread(datashade(xi_pz)) + dynspread(datashade(pr_pz))).cols(2)


