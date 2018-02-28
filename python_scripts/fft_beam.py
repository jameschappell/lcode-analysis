import h5py
import numpy as np
import matplotlib.pyplot as plt
import xi_r_map
from scipy import signal
import peakutils

n = 7.0e14
c = 2.9979245800e8
w_p = 2 * np.pi * 8980 * np.sqrt(n)
E_0 = 0.511e-3 / (c/w_p)
e = 1.6e-19

xi_min=-500.0
xi_max=-0.0
r_max=0.5

font = {'weight': 'normal', 'size': 15}
plt.rc('font', **font)

h5 = '/Users/jamiechappell/Documents/PhD/Simulations/lcode/multi_beam/' \
     'two_beam/gap_100/tb20000.h5'

Q, xi, x, y, px, py, pz = xi_r_map.h5import(h5)
r = np.sqrt(x*x + y*y)  # mm
H, xedges, yedges, xi_bin = xi_r_map.plot_xi_r_map(xi, r, Q, h5)
print H.shape


def FFT_analysis_beam_modulation(data):

    freq_array = []
    gamma_array = []

    for i in range(0, 1):

        width = H.shape[0]
        print width
        data = H[0][0:2499]
        #data = np.mean(H, axis=0)
        data_offaxis = H[0][2500:4999]

        fft_sample = xi_bin*1e-3 # m
        fft_sample = 1./(fft_sample / c) # s^-1
        n_len = len(data)

        data_xi = np.linspace(xi_min, xi_max, len(data))

        freqvec = 2 * np.pi * np.linspace(0, fft_sample/2, n_len/2)

        Y = np.fft.fft(data)/n
        Y = Y[range(n_len/2)]

        # off axis
        Z = np.fft.fft(data_offaxis)/n
        Z = Z[range(n_len/2)]

        indices = peakutils.indexes(abs(Y), thres=1e-12)
        indices_offaxis = peakutils.indexes(abs(Z), thres=1e-12)

        max_val = max(abs(Y)[indices])
        index = np.where(abs(Y) == max_val)
        pos_max = freqvec[index]

        # off axis
        max_val_oa = max(abs(Z)[indices_offaxis])
        index_oa = np.where(abs(Z) == max_val_oa)
        pos_max_oa = freqvec[index_oa]


        pos_string = '%s' % float('%.4g' % pos_max)
        str_max = '%s' % float('%.4g' % max_val)
        max_string = '($\omega$ = %s)' % pos_string

        # off-axis
        pos_string_oa = '%s' % float('%.4g' % pos_max_oa)
        max_string_oa = '($\omega$ = %s)' % pos_string_oa

        delta_omega = abs(pos_max - pos_max_oa)
        delta_string = '($\delta \omega$ = %.4g)' % delta_omega

        lambda_p = 1e3 * (2 * np.pi * c / pos_max) # mm
        gamma = ((pos_max**2 * 9.11e-31 * 8.85e-12)/((n*1e6) * e**2))**(-1)

        if gamma < 1:
            gamma = 1

        print 'Largest peak found at plasma frequency = %s s^-1' % pos_string
        print 'This corresponds to a plasma wavelength = %.4g mm' % lambda_p
        print 'An electron oscillating at this frequency has gamma = %.4g' % \
              gamma

        freq_array.append(pos_max)
        gamma_array.append(gamma)

        fig, ax = plt.subplots(3, 1)

        ax[0].pcolormesh(xedges, yedges, H, cmap='Blues', vmin=0, vmax=3e-12)
        ax[0].pcolormesh(xedges, -yedges, H, cmap='Blues', vmin=0, vmax=3e-12)
        ax[0].set_xlabel(r'$\xi\ (\mathrm{mm})$')
        ax[0].set_ylabel(r'$x\ (\mathrm{mm})$')
        ax[0].set_ylim(-r_max, r_max)
        ax[0].set_xlim(xi_min, xi_max)
        ax[1].plot(data_xi, data, 'r')
        ax[1].plot(data_xi, data_offaxis, 'g')
        ax[1].set_xlabel(r'$\xi$ $/mm$')
        #ax[1].set_ylabel(r'$N_p$')
        ax[1].set_xlim(xi_min, xi_max/2)
        ax[2].plot(freqvec, abs(Y), 'r', label='Second beam FFT')
        ax[2].plot(freqvec, abs(Z),'g', label='First beam FFT')
        ax[2].plot(pos_max, max_val, 'ro')
        ax[2].plot(pos_max_oa, max_val_oa, 'go')
        ax[2].annotate(delta_string, xy=(pos_max, 1.1 * max_val))
        ax[2].set_xlabel(r'$\omega$ $/s^{-1}$')#
        plt.legend()
        plt.show()

    return freq_array, gamma_array


freq, gamma = FFT_analysis_beam_modulation(H)

#plt.scatter(range(0, len(gamma)), gamma)
#plt.show()




