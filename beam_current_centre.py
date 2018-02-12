import matplotlib.pyplot as plt
import numpy as np
import os
import beamfile_output

n = 7.0e14
c = 2.9979245800e8
w_p = 2 * np.pi * 8980 * np.sqrt(n)
E_0 = 0.511e-3 / (c/w_p)

beamfileloc = '../awake-baseline/beamfile.bin'

h5 = os.path.splitext(beamfileloc)[0] + '.h5'

df0 = beamfile_output.read_beamfile(beamfileloc)
df = beamfile_output.calculate_xy(df0)



# selecting only protons
dfp = df[df.q > 0]

# setting r-limit
dfp = dfp[dfp.r_mm < 1.5]

dxi = 0.01 * (c/w_p)        # m (simulation step code used in lcode.cfg)
re = 2.818e-15              # m
q0 = 1.6e-19                # C
Q = ((dxi/re) * q0 / 2) * df['q'].head(1).values    # macroparticle charge Q

xi_min = -80.0
xi_max = 0.0
bin_width = 0.1         # mm

hist, bins = np.histogram(dfp['xi_mm'], range=(xi_min, xi_max),
                          bins=int((xi_max - xi_min) / bin_width))

I = hist * Q / (bin_width * 1e-3 / c)   # A

xi_centers = (bins[:-1] + bins[1:]) / 2

fig = plt.figure(figsize=(12, 4))

plt.bar(xi_centers, I, align='center', width=bin_width, edgecolor='none')

plt.xlim(xi_min, xi_max)
plt.ylim(0, 80)
plt.grid()
plt.xlabel(r'$\xi\ (\mathrm{mm})$')
plt.ylabel(r'$I\ (\mathrm{A})$')

plt.show()

