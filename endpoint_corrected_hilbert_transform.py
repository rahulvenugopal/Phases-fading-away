#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 11:07:35 2022

@author: Sruthi Kuriakose

The original code in Matlab - Supplementary code of 
Schreglmann, Sebastian R et al. “Non-invasive suppression of essential tremor 
via phase-locked disruption of its temporal coherence.” Nature communications 
vol. 12,1 363. 13 Jan. 2021, doi:10.1038/s41467-020-20581-7

"""
import numpy as np
import scipy.signal
from scipy import fft as sp_fft
import matplotlib.pyplot as plt


#%%

def echt_scipy_hilbert(x, freq,samplingrate, N=None, axis=-1):
    '''scipy implementation of hilbert (scipy/signal/_signaltools.py) + echt'''
    '''
    Parameters
    ----------
    x : array_like
        Signal data.  Must be real.
    freq : int
        Frequency of signal - used to apply causal filter
    samplingrate : int
        Sampling rate of signal
    N : int, optional
        Number of Fourier components.  Default: ``x.shape[axis]``
    axis : int, optional
        Axis along which to do the transformation.  Default: -1.

    Returns
    -------
    x : ndarray
        Analytic signal of `x`, of each 1-D array along `axis`
    '''
    x = np.asarray(x)
    nsamples = len(x)
    if np.iscomplexobj(x):
        raise ValueError("x must be real.")
    if N is None:
        N = x.shape[axis]
    if N <= 0:
        raise ValueError("N must be positive.")

    Xf = sp_fft.fft(x, N, axis=axis)
    h = np.zeros(N, dtype=Xf.dtype)
    if N % 2 == 0:
        h[0] = h[N // 2] = 1
        h[1:N // 2] = 2
    else:
        h[0] = 1
        h[1:(N + 1) // 2] = 2

    if x.ndim > 1:
        ind = [np.newaxis] * x.ndim
        ind[axis] = slice(None)
        h = h[tuple(ind)]
    Xf = Xf*h
    
    # 'echt'
    # Compute filter's frequency response
    filt_lf = freq - freq/4
    filt_hf = freq + freq/4
    filt_order = 2;
    b, a = scipy.signal.butter(filt_order, [filt_lf,filt_hf], 'bandpass', analog=False,fs=samplingrate)
    # Center bin frequencies of the FFT.
    filt_freq = np.ceil(np.arange(-nsamples/2, nsamples/2))
    #Freq domain filter coefficents.
    filt_coeff = scipy.signal.freqz(b,a,worN = filt_freq,fs = samplingrate)[1]

    
    fftdata = Xf
    # % Multiply FFT by filter's response function 
    fftdata = np.fft.fftshift(fftdata)  # rearranges array so DC at the center
    fftdata = fftdata * filt_coeff;  
    fftdata = np.fft.ifftshift(fftdata);   # rearranges back array so DC at the left
    
    #IFFT
    x = sp_fft.ifft(fftdata, axis=axis)
    return x


def echt_rdgao_hilbert(x,freq,samplingrate):  
    ''' rdgao's implementation of hilbert from http://www.rdgao.com/roemerhasit_Hilbert_Transform/ 
        + echt '''
    '''
    Parameters
    ----------
    x : array_like
        Signal data.  Must be real.
    freq : int
        Frequency of signal - used to apply causal filter
    samplingrate : int
        Sampling rate of signal
    Returns
    -------
    x : ndarray
        Analytic signal of `x`, of each 1-D array along `axis`
    '''
    nsamples = len(x)
    Fx = np.fft.fft(x) # get fft
    f_axis = np.fft.fftfreq(len(x))
    Fx[np.where(f_axis<0)]=0. # zero out negative frequencies
    Fx=Fx*2 # return 2x positive frequencies
    
    
    # % Compute filter's frequency response
    filt_lf = freq - freq/4
    filt_hf = freq + freq/4
    filt_order = 2;
    b, a = scipy.signal.butter(filt_order, [filt_lf,filt_hf], 'bandpass', analog=False,fs=samplingrate)
    # Center bin frequencies of the FFT.
    filt_freq = np.arange(-nsamples/2, nsamples/2)
    #Freq domain filter coefficents.
    filt_coeff = scipy.signal.freqz(b,a,worN = filt_freq,fs = samplingrate)[1]

    
    fftdata = Fx
    # % Multiply FFT by filter's response function 
    fftdata = np.fft.fftshift(fftdata)  # rearranges array so DC at the center
    
    fftdata = fftdata * filt_coeff;  

    fftdata = np.fft.ifftshift(fftdata);   # rearranges back array so DC at the left
   
    # IFFT
    z= np.fft.ifft(fftdata) 
    
    return z



