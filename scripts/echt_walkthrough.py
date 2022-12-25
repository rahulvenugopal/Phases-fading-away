# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 07:35:40 2022
- Create pure sine waves of specific frequency
- Create many wave forms where their starting and ending phases vary
- Replicate the phase distortions due to Gibbs phenomenon
- Cure those distortions by end point corrected Hilbert magic
@author: Rahul Venugopal
"""
#%% Libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft
from scipy import signal

# Sine wave generation
sampling_rate = 1000

# creating sine wave time series where the starting phase is different from ending phase
# this is where FFT and Hilbert inherits phase errors
length_seconds = [2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4]
fi = 2

# Different phases where wave starts, have used slightly longer length_seconds
# so that start and end phases are different

phase = [0, np.pi/8, np.pi/4, np.pi/2, np.pi, -np.pi/8, -np.pi/4, -np.pi/2, -np.pi]

for lengths in length_seconds:

    for entries in phase:

        npnts = sampling_rate*lengths  # number of time samples

        time = np.arange(0, npnts)/sampling_rate

        signalu = np.around(np.sin(2*np.pi*fi* (time - entries)), decimals = 2)

        # mins = [i for i, x in enumerate(signalu) if x == min(signalu)]

        # plotting
        plt.figure(figsize=(15,10))
        plt.subplot(6,1,1)
        plt.plot(time,signalu, 'k', label='Raw Signal', alpha=0.4, lw=1)

        plt.subplot(6,1,2)
        fft_signal = fft(signalu)
        plt.plot(fft_signal)

        # Analytic signal and phase instantaneous
        analytic_signal = signal.hilbert(signalu)
        phase_diff = np.angle(analytic_signal)

        # plotting
        plt.subplot(6,1,3)
        plt.plot(phase_diff, 'indianred', ms=1)
        plt.yticks([])
        plt.ylabel('Instantaneous phase')

        # add vertical lines | Commenting it out now, gotta do it neatly
        # [plt.axvline(_x, color="#d73027", ls = 'dashdot', linewidth = 0.1) for _x in mins]

        # plotting real part of analytical signal
        plt.subplot(6,1,4)
        plt.plot(analytic_signal.real, 'steelblue', ms=1)
        plt.yticks([])
        plt.ylabel('Real part')

        # plotting imaginary part of analytical signal
        plt.subplot(6,1,5)
        plt.plot(analytic_signal.imag, 'steelblue', ms=1)
        plt.yticks([])
        plt.ylabel('Imag part')

        # plotting envelope
        plt.subplot(6,1,6)
        plt.plot(np.abs(analytic_signal), 'steelblue', ms=1)
        plt.yticks([])
        plt.ylabel('Envelope')

        plt.tight_layout()

        plt.savefig(str(lengths) + '_seconds_' + str(entries) + '_Hilbert phases.png',
                    dpi = 600,
                    backend='cairo')
        plt.close()
        print(str(entries) + '_' + str(lengths))