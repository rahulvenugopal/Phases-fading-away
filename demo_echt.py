# -*- coding: utf-8 -*-
"""
Created on Sat Dec 24 15:23:01 2022

@author: Sruthi Kuriakose
"""
#%% Load libraries
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt

# Load the script with custom functions
# You need to be in the correct path to do this
import echt as echt

#%% replica of signal from paper
samplingrate = Fs = 2048  # sampling freq
freq = 2.25 #2.25 Hz frequency
T = 1#seconds

nsamples = Fs* T
t = np.arange(nsamples)/Fs

signal = np.cos(2 * np.pi * freq * t)
plt.plot(t, signal,label='true_signal')

# true_imag_signal
true_imag_signal = np.cos(2 * np.pi * freq * t - np.pi/2)

plt.plot(t,true_imag_signal,label='signal phase-shifted by pi/2')

# converts to complex signal
scipy_analytic_signal = scipy.signal.hilbert(signal)

scipy_hilbert_data = scipy_analytic_signal.imag
plt.plot(t,scipy_hilbert_data, label='scipy hilbert')
plt.legend()

#%% Let the endpoint correction happen
xsignal = signal.copy()

# z = echt.echt_rdgao_hilbert(xsignal, freq)
# plt.plot(t,z.imag,label='rdgao + echt')

z = echt.echt_scipy_hilbert(xsignal, freq,samplingrate)
plt.plot(t,z.imag,label='scipy hilbert + echt')

plt.legend()