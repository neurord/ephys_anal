import numpy as np
from matplotlib import pyplot as plt

def next2pow(x):
    return 2**int(np.ceil(np.log(float(x))/np.log(2.0)))

def rolling_window(data,window_size,step):
    num_steps=1+(len(data)-window_size)/step
    begin=[n*step for n in range(num_steps)]
    slice=np.zeros((num_steps,window_size))
    for i,beg in enumerate(begin):
        slice[i,:]=data[beg:beg+window_size]
    return slice
    
def fft_func(wave,ts):
    fft_wave = fftpack.fft(wave)
    #Note that maximum frequency for fft is fs/2; the frequency unit is cycles/time units.
    #freqs is x axis. Two ways to obtain correct frequencies:
    #specify sampling spacing as 2nd parameter, note that fs=1/ts
    #freqs = fftpack.fftfreq(len(wave))*fs
    #multiply by sample spacing by max frequency (=1/ts):
    freqs=fftpack.fftfreq(len(wave),ts)
    phase=np.arctan2(fft_wave.imag,fft_wave.real)
    return fft_wave,phase,freqs

plt.ion()
def fft_plot(t,x,freqs,X,phase=None,title=None):
    if np.array(phase).any():
        fig,axes=plt.subplots(3,1)
    else:
        fig,axes=plt.subplots(2,1)
    fig.suptitle(title)
    fs=np.abs(np.min(freqs))
    axes[0].plot(t,x)
    axes[0].set_xlabel('time')
    axes[0].set_xlim(0,t[-1])

    axes[1].plot(freqs, np.abs(X),'.')
    axes[1].set_xlabel('Frequency in Hertz [Hz]')
    axes[1].set_ylabel('FFT Magnitude')
    axes[1].set_xlim(-fs , fs )
    
    if np.array(phase).any():
        axes[2].plot(freqs,phase,'r.')
        axes[2].set_xlabel('Frequency in Hertz [Hz]')
        axes[2].set_ylabel('FFT Phase')
        axes[2].set_xlim(-fs , fs )
    return

def plot_spectrum(spectrum,phase_spect,S_db):
    ############## No color bar in this version
    f,axes = plt.subplots(nrows=3,ncols=1,figsize=(4.8, 4.8))
    f.suptitle('spectrum')
    limits=[0, tend, 0, fs / 2 ]
    im0=axes[0].imshow(S_db, origin='lower', cmap='viridis',
                       extent=limits)
    im1=axes[1].imshow(np.abs(spectrum), origin='lower', cmap='viridis',
                       extent=limits)
    im2=axes[2].imshow(phase_spect, origin='lower', cmap='viridis',
                       extent=limits)
    #f.colorbar(im1,cax)
    for ax in axes:
        ax.axis('tight')
        ax.set_ylabel('Frequency [kHz]')

    ax.set_xlabel('Time [s]')

    ############## Color bar works in this version, which is uglier
    plt.figure()
    plt.suptitle('spectrum')

    plt.subplot(311)
    plt.imshow(S_db, origin='lower', cmap='viridis', aspect='auto',extent=limits)
    plt.colorbar()
    plt.ylabel('Frequency [kHz]')

    plt.subplot(312)
    plt.imshow(np.abs(spectrum), origin='lower', cmap='viridis', aspect='auto',extent=limits)
    plt.colorbar()
    plt.ylabel('Frequency [kHz]')

    plt.subplot(313)
    plt.imshow(phase_spect, origin='lower', cmap='viridis', aspect='auto',extent=limits)
    plt.colorbar()
    plt.ylabel('Frequency [kHz]')

    plt.xlabel('Time [s]')
