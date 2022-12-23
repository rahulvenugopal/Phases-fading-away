function x = echt(xr, filt_lf, filt_hf, Fs, n)
% ECHT: Endpoint Corrected Hilbert Transform
%   X = echt(Xr) computes the so-called discrete-time analytic signal. 
%   It is a modification of Matlab's hilbert.m function that significantly 
%   ameliorates the distortions, known as Gibbs phenomenon, at the end of 
%   the signal.
%
%   NOTES:
%       1. One common implementation of the Hilbert Transform uses a DFT
%       (aka FFT) as part of its computation. Inherent to the DFT is the 
%       assumption that a finite sample of a signal is replicated 
%       infinitely in time, effectively abutting the end of a sample with
%       its replicated start. If the start and end of the sample are not
%       continuous with each other, distortions are introduced by the DFT.
%       Echt effectively smooths out this 'discontinuity' by 
%       selectively deforming the start of the sample. It is hence most 
%       suited for real-time applications in which the point/s of interest 
%       is/are the most recent one/s (i.e. last) in the sample window. 
%       2. We found that a filter bandwidth (BW=filt_hf-filt_lf) of up to
%       half the signal's central frequency works well.
%   
%   INPUT PARAMETERS:
%       xr: time domain signal 
%       filt_lf: low-cutoff frequency of a bandpass causal filter
%       filt_hf: high-cutoff frequency of a bandpass causal filter
%       Fs: sampling rate of time domain signal 
%       n: length of analytic signal (optional)
%
%   USE EXAMPLE:
%       f0 = 2;
%       filt_BW = f0/2;
%       N = 1000;
%       Fs = N/(2*pi);
%       t = -2*pi:1/Fs:0;
%       Xr = cos(2*pi*f0*t-pi/4);
%       Filt_LF = f0 - filt_BW/2;
%       Filt_HF = f0 + filt_BW/2;
%       x = echt(Xr, Filt_LF, Filt_HF, Fs, length(Xr))
%
%   See also: hilbert, FFT, IFFT, butter, freqz
%
%   References:
%       S. R. Schreglmann1*, D. Wang*, R. Peach*, J. Li, X. Zhang, E. Panella, 
%       E. S. Boyden, M. Barahona, S. Santaniello, K. P. Bhatia, J. Rothwel, N. Grossman
%       "Non-invasive Amelioration of Essential Tremor via Phase-Locked
%       Disruption of its Temporal Coherence".
%       * These authors contributed equally to this work

% Check input
if nargin<5, n=[]; end
if ~isreal(xr)
  warning(message('signal:hilbert:Ignore'))
  xr = real(xr);
end

%  Work along the first nonsingleton dimension
[xr,nshifts] = shiftdim(xr); 
if isempty(n)
  n = size(xr,1);
end

% Compute FFT
x = fft(xr,n,1); % n-point FFT over columns.

% Set negative components to zero and multiply by 2 postive (apart from DC
% and Nyquist frequency)
h  = zeros(n,~isempty(x)); % nx1 for nonempty. 0x0 for empty.
if n > 0 && 2*fix(n/2) == n
  % even and nonempty
  h([1 n/2+1]) = 1;
  h(2:n/2) = 2;
elseif n>0
  % odd and nonempty
  h(1) = 1;
  h(2:(n+1)/2) = 2;
end
x = x.*h(:,ones(1,size(x,2)));

% Compute filter's frequency response
filt_order = 2;
[b,a] = butter(filt_order, [filt_lf filt_hf]/(Fs/2), 'bandpass');
T = 1/Fs*n;                           % Sampling time period
filt_freq = (ceil(-n/2:n/2-1)'/T);    % Center bin frequencies of the FFT.
filt_coeff = freqz(b,a,filt_freq,Fs); % Freq domain filter coefficents.

% Multiply FFT by filter's response function 
x = fftshift(x);    % rearranges array so DC at the center
x = x.*filt_coeff;  
x = ifftshift(x);   % rearranges back array so DC at the left

% IFFT
x = ifft(x);

% Convert back to the original shape.
x = shiftdim(x,-nshifts);
