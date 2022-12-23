# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 07:35:40 2022

@author: Rahul Venugopal
"""
#%% Libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
from scipy import signal

# Sine wave generation
sampling_rate = 1000

# creating sine wave time series where the starting phase is different from ending phase
# this is where FFT and Hilbert inherits phase errors
length_seconds = [2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4]
fi = 2


phase = [0, np.pi/8, np.pi/4, np.pi/2, np.pi, -np.pi/8, -np.pi/4, -np.pi/2, -np.pi]

for lengths in length_seconds:
    
    for entries in phase:   
        
        npnts = sampling_rate*lengths  # number of time samples
    
        time = np.arange(0, npnts)/sampling_rate
        
        signalu = np.around(np.sin(2*np.pi*fi* (time - entries)), decimals = 2)
        
        mins = [i for i, x in enumerate(signalu) if x == min(signalu)]
        
        # plotting
        plt.figure(figsize=(15,6))
        plt.subplot(3,1,1)
        plt.plot(time,signalu, 'k', label='Raw Signal', alpha=0.4, lw=1)
        plt.ylabel('Real Signal')
        
        plt.subplot(3,1,2)
        fft_signal = fft(signalu)
        plt.plot(fft_signal)
        
        phase_diff = np.angle(signal.hilbert(signalu))
    
        # plotting
        plt.subplot(3,1,3)
        plt.plot(phase_diff, '.k', ms=1)
        plt.yticks([])
        plt.ylabel('Instantaneous phase')
        
        # add vertical lines
        [plt.axvline(_x, color="#d73027", ls = 'dashdot', linewidth = 0.05) for _x in mins]
    
        plt.tight_layout()
        
        plt.savefig(str(lengths) + '_seconds_' + str(entries) + '_Hilbert phases.png',
                    dpi = 600,
                    backend='cairo')
        plt.close()
        print(str(entries) + '_' + str(lengths))
