# -*- coding: utf-8 -*-
"""
Created on Sat Dec 24 15:23:01 2022

@author: Sruthi Kuriakose
"""
import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
import endpoint_corrected_hilbert_transform as echt

# replica of signal from paper
samplingrate = Fs = 2047  # sampling freq
freq = 2.25 #2.25 #frequency
T = 1#seconds
# T=1/samplingrate*nsamples  #sampling time period

nsamples = Fs* T
t = np.arange(nsamples)/Fs
signal = np.cos(2 * np.pi * freq * t)

plt.plot(t, signal,label='true_signal')

# true_imag_signal = np.sin(2 * np.pi * freq * t)
true_imag_signal = np.cos(2 * np.pi * freq * t - np.pi/2)
# when phase = np.pi/2, and siganl was np.sin(2 * np.pi * f * t + phase), true_imag_signal = np.sin((2 * np.pi * f * t + phase) - np.pi/2  )

plt.plot(t,true_imag_signal,label='signal phase-shifted by pi/2')

scipy_analytic_signal = scipy.signal.hilbert(signal)  # converts to complex signal
scipy_hilbert_data = scipy_analytic_signal.imag
plt.plot(t,scipy_hilbert_data, label='scipy hilbert')

plt.legend()

#%%
xsignal = signal.copy()

# z = echt.echt_rdgao_hilbert(xsignal, freq)
# plt.plot(t,z.imag,label='rdgao + echt')

z = echt.echt_scipy_hilbert(xsignal, freq,samplingrate)
plt.plot(t,z.imag,label='scipy hilbert + echt')

plt.legend()

#%%