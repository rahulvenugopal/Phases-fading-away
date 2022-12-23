# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 17:24:25 2022
- Understand the phase distortion explained to us (Sruthi, Aiswarya, Rahul) by
Nir Grossman at EMBO lecture course IIT Gn 2022

- The whole idea was that Hilbert transform (the way it is done) introduces an
error in the endpoint phase, this is in the range of 160 degrees or so (worst case)
- The error goes away when the start and end are mimicking each other
@author: Rahul Venugopal
"""
#%% Richard Gao tutorial
# http://www.rdgao.com/roemerhasit_Hilbert_Transform/

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

plt.rcParams["axes.labelsize"]=20
plt.rcParams["font.size"]=15
plt.rcParams["font.family"] = "Arial"
from scipy import io, signal

# our open source lab package, free for all to use/contribute :)
# https://github.com/voytekresearch/neurodsp
import neurodsp as ndsp
from neurodsp.filt import filter_signal

#%%

# loading the test data
data = np.load('sample_data_2.npy')
x = data[:20000]
fs = 1000.
t = np.arange(len(x))/fs

# filtering and hilbert transform
x_filt = filter_signal(x, fs, 'bandpass', f_range=(4,12), remove_edges=False)
x_a = signal.hilbert(x_filt)

# plotting
plt.figure(figsize=(15,6))
plt.subplot(2,1,1)
plt.plot(t,x, 'k', label='Raw Signal', alpha=0.4, lw=1)
plt.plot(t,x_filt, 'r', label='Filtered Signal', lw=2)
plt.ylabel('Real Signal')
plt.legend()
plt.xlim((0,3));

plt.subplot(2,1,2)
plt.plot(t,x_filt, '--r', label='Filtered Signal', alpha=1, lw=1)
plt.plot(t,x_a.real, label='Hilbert Real', alpha=0.5, lw=2)
plt.plot(t,x_a.imag, label='Hilbert Imag', alpha=0.5, lw=2)
plt.xlabel('Time (s)')
plt.ylabel('Complex Signal')
plt.legend(loc='upper left')
plt.xlim((0,3));

# selecting 2 slices of data
t1 = 2650
t2 = 2700
yb = plt.ylim()
plt.plot([t1/fs]*2, yb, 'k--')
plt.plot([t2/fs]*2, yb, 'g--')
plt.tight_layout()

#%%
x_power = np.abs(x_a)**2
x_phase = np.angle(x_a)
# plotting
plt.figure(figsize=(15,6))
plt.subplot(2,1,1)
plt.plot(t,x_power, '.k', ms=1)
plt.yticks([])
plt.ylabel('Power')
plt.xlim((0,3));

plt.subplot(2,1,2)
plt.plot(t,x_phase, '.k', ms=1)
plt.ylabel('Phase')
plt.xlabel('Time (s)')
yb = plt.ylim()
plt.plot([t1/fs]*2, yb, 'k--')
plt.plot([t2/fs]*2, yb, 'g--')
plt.xlim((0,3))
plt.tight_layout()

#%%
plt.figure(figsize=(6,6))
plt.axhline(color='k', lw=1)
plt.axvline(color='k', lw=1)
plt.plot(x_a.real[:5000],x_a.imag[:5000], 'k.-', alpha=0.1, ms=2)
plt.quiver(0,0,x_a[t1].real,x_a[t1].imag, angles='xy', scale_units='xy', scale=1, color='k')
plt.quiver(0,0,x_a[t2].real,x_a[t2].imag, angles='xy', scale_units='xy', scale=1, color='g')
plt.xlabel('Real');plt.ylabel('Imaginary')
plt.box('off')
plt.title('Vector 1(k) phase: %.2f, vector 2(g) phase: %.2f'%(x_phase[t1]/(2*np.pi)*360, x_phase[t2]/(2*np.pi)*360));

#%% Synthetic signals
# create a sine wave for 5 seconds at 14Hz, with slight phase delay
t_trunc = t[:5000]
x_sin = np.cos(2*np.pi*14*t_trunc+0.5)
F_sin = np.fft.fft(x_sin, norm='ortho')
f_axis = np.fft.fftfreq(len(x_sin),1/fs)

plt.figure(figsize=(15,10))
# time domain plot
plt.subplot(4,1,1)
plt.plot(t_trunc,x_sin)
plt.xlabel('Time (s)')
plt.ylabel('Real Signal')
plt.xlim((0,3))


# frequency domain plots
plt.subplot(4,1,2)
plt.plot(f_axis, np.abs(F_sin)**2,'ko-')
plt.ylabel('Power')
plt.title('Fourier Power')
plt.xlim([-500,500])
plt.subplot(4,1,3)
plt.plot(f_axis, F_sin.real, 'ko-')
plt.title('Fourier Coeff. Real')
plt.xlim([-500,500])
plt.subplot(4,1,4)
plt.plot(f_axis, F_sin.imag, 'ko-', label='Fourier Imag')
plt.title('Fourier Coeff. Imag')
plt.xlim([-500,500])
plt.xlabel('Frequency (Hz)')
plt.tight_layout()

#%% Create test sine waves

# Libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.fft import fft, fftfreq

phase = [0, np.pi/8, np.pi/4, np.pi/2, np.pi, -np.pi/8, -np.pi/4, -np.pi/2, -np.pi]

for entries in phase:

    time = np.arange(12000)*np.pi/1000
    signalu = np.sin(time - entries)

    # FFT
    signalu_fft = fft(signalu)

    xf = fftfreq(12000, 1/1000)[:12000//2]
    plt.plot(xf, 2.0/12000 * np.abs(signalu_fft[0:12000//2]))

    # Hilbert
    x_a = signal.hilbert(signalu)

    # plotting
    plt.figure(figsize=(15,6))
    plt.subplot(3,1,1)
    plt.plot(time,signalu, 'k', label='Raw Signal', alpha=0.4, lw=1)
    plt.ylabel('Real Signal')
    plt.legend()


    plt.subplot(3,1,2)
    plt.plot(time,signalu, '--r', label='Filtered Signal', alpha=1, lw=1)
    plt.plot(time,x_a.real, label='Hilbert Real', alpha=0.5, lw=2)
    plt.plot(time,x_a.imag, label='Hilbert Imag', alpha=0.5, lw=2)
    plt.xlabel('Time (s)')
    plt.ylabel('Complex Signal')
    plt.legend(loc='upper left')
    plt.show()
    
    phase_diff = np.angle(x_a)

    # plotting
    plt.subplot(3,1,3)
    plt.plot(phase_diff, '.k', ms=1)
    plt.yticks([])
    plt.ylabel('Phase difference between real and imag')

    plt.tight_layout()

    plt.savefig(str(entries) + 'Hilbert phases.png',
                dpi = 600,
                backend='cairo')
    plt.close()

# Angle bussiness
phase_diff = np.angle(x_a)

# plotting
plt.figure(figsize=(15,6))
plt.plot(phase_diff, '.k', ms=1)
plt.yticks([])
plt.ylabel('Phase difference between real and imag')

plt.tight_layout()