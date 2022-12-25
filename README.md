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
- We created a sine wave of `2.5 Hz` with a sampling frequency of `1000 Hz`. $${\color{blue} Blue curve}$$
- A phase shifted waveform was generated which serves as the ground truth. $${\color{orange} Orange curve}$$
- Hilbert transform was run on the sine wave and we took out the imaginary part. $${\color{green} Green curve}$$
- We applied the endpoint corrected Hilbert and the distortions to phase (*stemming from distortions to the imaginary part of analytical signal*) got corrected. $${\color{red} Red curve}$$
- Look at the price we paid for this, the beginning of the red curve suffers from a phase shift!

$${\color{green} Green curve}$$

![](https://github.com/rahulvenugopal/Phases-fading-away/blob/main/Demo.png)


