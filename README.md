# Phases fading away
We (Sruthi, Aiswarya and Rahul) met Dr. Nir Grossman at EMBO Lecture course on `Neuromodulation` at IIT Gandhinagar and he walked us through endpoint Hilbert method under :flashlight: (because we were explaining our poster around 7 PM and Dr. Grossman was really patient and gave us valubale feedbacks) Dr. Nir, you are awesome :fire:

## This repo aims to understand what he taught us in detail
- If the start and end of the sample are not continuous with each other, distortions are introduced by the DFT
- Echt effectively smooths out this 'discontinuity' by selectively deforming the start of the sample (**there is no free lunch**)
- We wanted to capture realtime data from EEG, run filter-Hilbert to find the instantaneous phase and trigger sound if the phase is of certain value
- Since the conventional method distorts phases at the end point, we need to fix this
- The redboxes are distortions! Damn
- Gibbs Gibbs Gibbs :bomb:

## Demo
- `demo_echt.py` contains the demo script and the Pythonified code for endpoint correction can be found under scripts folder `echt.py`
- We created a sine wave of `2.5 Hz` with a sampling frequency of `1000 Hz`. Blue curve
- A phase shifted waveform was generated which serves as the ground truth. Orange curve
- Hilbert transform was run on the sine wave and we took out the imaginary part. Green curve
- We applied the endpoint corrected Hilbert and the distortions to phase (*stemming from distortions to the imaginary part of analytical signal*) got corrected. Red curve
- Look at the price we paid for this, the beginning of the red curve suffers from a phase shift!
![](https://github.com/rahulvenugopal/Phases-fading-away/blob/main/Demo.png)

- Let us understand, how endpoint correction really works

## To Do
1. Replicate the phase and amplitude error plots from [Schreglmann, S.R et al 2021](https://www.nature.com/articles/s41467-020-20581-7#citeas)
- The maximum phase error was for 90 degrees and 270 degress the error was minimal!!!! Understand why
2. Make sure the echt.py implementation is correct by bench-marking against original MATLAB implementation
- We spotted the `filt_freq` need to be adjusted based on the length of the samples (even or odd) using a `np.ceil` function
3. Understand what causal filter is doing to the FFT in an intuitive sense
- From the `Bressler et al 2022` paper:
> the amplitude and phase at the beginning of the sample window is warped to be continuous with the end

4. Read [A Wearable EEG System for Closed-Loop Neuromodulation of High-Frequency Sleep-Related Oscillations](https://arxiv.org/abs/2212.11273)
5. Create neat visualisation for phase error plots as a circular histogram with custom annotations
![Image from Figure 7 of Bressler et al, 2022 paper](https://github.com/rahulvenugopal/Phases-fading-away/blob/main/Circular_barhistogram_with_custom_texts.png)

## Explorations
1. Understand the causal filter design in an intuitive manner
2. Delay the real time fetching so that the target portion does not fall on distortion part and target the next onset of phase
