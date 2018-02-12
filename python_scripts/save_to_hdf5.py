import h5py
import beamfile_output
import numpy as np
import os

n = 7.0e14
c = 2.9979245800e8
w_p = 2 * np.pi * 8980 * np.sqrt(n)
E_0 = 0.511e-3 / (c/w_p)

numbers = np.arange(34000, 50001, 2000)
file_list = []
location = '/unix/pdpwa/jchappell/lcode/multi_beam/two_beam' \
              '/tb'

for i in numbers:
    number_string = '{0:05}'.format(i)
    string = location + number_string + '.swp'
    file_list.append(string)


for beamfileloc in file_list:

    h5 = os.path.splitext(beamfileloc)[0] + '.h5'

    df0 = beamfile_output.read_beamfile(beamfileloc)
    df = beamfile_output.calculate_xy(df0)

    dxi = 0.01 * (c/w_p)        # m (simulation step code used in lcode.cfg)
    re = 2.818e-15              # m
    q0 = 1.6e-19                # C
    Q = ((dxi/re) * q0 / 2) * df['q'].head(1).values    # macroparticle charge Q

    print("Saving beamfile to " + h5 + "...")
    os.path.exists(h5) and os.remove(h5)

    with h5py.File(h5, 'w') as hf:
        hf.create_dataset("macroparticle_charge", data=Q, dtype='float32') # in C

        for col in ['x_mm', 'y_mm', 'xi_mm', 'px_GeV', 'py_GeV', 'pz_GeV']:

            print(col+' ')
            ds = hf.create_dataset(col, data=df[col].values, dtype='float32')

print("\nDone.")