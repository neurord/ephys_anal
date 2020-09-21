import numpy as np
from scipy import fftpack, signal
import matplotlib.pyplot as plt
import fft_utils as u

############################
#           parameters
############################
#20 s of signal
tend=20
#sampling interval
ts=0.01
#sampling frequency (measurements per second)
f0 = 1 #starting frequency, Esprit used 10
f1 = 10 #ending frequency, Esprit used 40
desired_fs=10*f1

#better sampling frequency - power of two closest to desired_fs
fs=u.next2pow(desired_fs)
ts=1./fs

t = np.linspace(0, tend, fs*tend)


#########################################
## create and analyze sum of sinusoids
#########################################

phi0=0
phi1=np.pi/2
A0=1
A1=2
sine0 = A0*np.sin(f0 * 2 * np.pi * t+phi0)
sine1 = A1*np.sin(f1 * 2 * np.pi * t+phi1)

fft0,phase0,freq0=u.fft_func(sine0,ts)
fft1,phase1,freq0=u.fft_func(sine1,ts)
fft_tot,phasetot,freq0=u.fft_func(sine0+sine1,ts)
#Note that phase is non-zero even where FFT magnitude is 0.  What does that mean?

#plot results
u.fft_plot(t,sine0,freq0,fft0,phase=phase0)
u.fft_plot(t,sine1,freq0,fft1,phase=phase1)
u.fft_plot(t,sine1+sine0,freq0,fft_tot,phase=phasetot)

#########################################
## create and analyze chirp using DFT
#########################################

chrp = signal.chirp(t, f0, tend, f1, method='linear')
#chirp = cos[2*pi(f0*t+(k/2)*t^2)
#k=(f1-f0)/T, where T is duration (tend)
fft_chrp,chrp_phase,freq1=u.fft_func(chrp,ts)

u.fft_plot(t,chrp,freq1,fft_chrp,phase=chrp_phase,title='DFT')

#########################################
#     Repeat analysis using chirpz transform
# https://github.com/ericmjonas/pychirpz/blob/master/chirpz/pychirpz.py
#########################################
import pychirpz as czt

m = len(chrp)                          ## number of points desired
w = np.exp(-2j*np.pi*(f1-f0)/(m*fs))  ## freq. step of f2-f1/m
a = np.exp(2j*np.pi*f0/fs)
czt_chirp=czt.chirpz(chrp,m,w,a)

phase_czt = np.arctan2(czt_chirp.imag, czt_chirp.real)
u.fft_plot(t,chrp,freq1,czt_chirp,phase=phase_czt,title='chirpz')

#########################################
#         spectrogram using rolling window
#https://www.oreilly.com/library/view/elegant-scipy/9781491922927/ch04.html
#########################################
M=512 #size of window - arbitrary
offset=100 #each slice is offset by 100 from the next slice
data=chrp
slices=u.rolling_window(data,M, offset) 
win = np.hanning(M + 1)[:-1] #filter the window to prevent edge or wrap-around effects
slices = slices * win
slices = slices.T
#one-sided spectrum from M // 2 + 1
spectrum = fftpack.fft(slices, axis=0)[:M // 2 + 1:-1]
phase_spect=np.arctan2(spectrum.imag, spectrum.real)
#convert to decibels
S_db = 20 * np.log10(np.abs(spectrum) / np.max(np.abs(spectrum)))
plot_spectrum(spectrum,phase_spect,S_db)

## Calculation of impedance
''' Need both input signal and output signal
calculate fft of both input and output using (A) DFT, (B) chirpz, or (C) spectrum
Real part of impedance (resistance) is np.abs(fft output)/np.abs(fft input)
Complex part of impedance is phase of output - phase of input
Eqn: Z=V0/I0 * e^i*phase
Re(Z) = np.abs(fftV)/np.abs(fftI)
Im(Z) = phase(fftV)-phase(fftI)
likely the same thing should be done for output of chirpz
'''


