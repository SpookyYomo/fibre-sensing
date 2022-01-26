# Initialiastion: Directory appending for my system. Vary the directories as necessary.
import sys, os.path

from numpy.core.fromnumeric import argmax
# Add the FYPLibrary folder in the same level as this folder.
if os.path.dirname(os.path.dirname(__file__)) not in sys.path:
    # print(f"{os.path.dirname(os.path.dirname(__file__)) = }")
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'FYPLibrary'))

from pprint import pprint
# pprint(sys.path)

# Import Modules
from os import listdir
from file_reading import *
from IQ_demod import *
import numpy as np
import EPstandard
from scipy.signal import periodogram
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# Yield file in batch folder
# folderpath = os.path.dirname(__file__)+'\Batch 2b'
# files = listdir(path=folderpath)
# NAME = files[5]
folderpath = os.path.dirname(__file__)+'\periodogram_test'
files = listdir(path=folderpath)
NAME = files[0]
TXT_FILE_PATH = file_str_to_path(NAME, dir=folderpath)
print(f"{TXT_FILE_PATH = }")
_, mirror_f, mVpp, delta_corr = numerical_variables_from_name(NAME)

SIGNAL_F = 80.125e6*2 #Hz 
SAMPLING_F = 1.0e6 #Hz
N, M = freq_ratio(signal=SIGNAL_F, sample=SAMPLING_F)

# file reading
meta, trace = fr.parse_and_read_oscilliscope_txt(TXT_FILE_PATH)
signal = signal_from_trace(np.asarray(trace))
phases = signal_to_phase(signal, N, 2*pi/N, True)
phases = phase_reconstruction(phases, 4.8)
amplitudes = get_R_signal(signal, N, 2*pi/N) # amplitude over time
t_axis = np.arange(start= 0, 
    stop= (int(meta["Record Length"][0])-N) * meta['Sample Interval'][0], step= meta['Sample Interval'][0])
time_axis = np.arange(start= 0, 
    stop= (int(meta["Record Length"][0])) * meta['Sample Interval'][0], step= meta['Sample Interval'][0])

# generate plot
fig, axs = plt.subplots(nrows = 2, ncols= 3, sharex= False)

# signal plot
ax = axs[0, 0]
ax.plot(time_axis, signal, color = 'mediumblue')
ax.set_ylabel(r'$y(t)$/V', useTex = True)
ax.set_xlabel(r'$t$/s', useTex = True)
ax.set_title(f"Signal", useTex = False)

# phase plot
ax = axs[0, 1]
ax.plot(t_axis, phases, color = 'mediumblue')
# ax.plot(t_axis, linear(t_axis, *popt), color = 'darkgreen', markersize = 3)
ax.set_ylabel(r'$\phi$/rad', useTex = True)
ax.set_xlabel(r'$t$/s', useTex = True)
ax.set_title(f"Phase", useTex = False)

# amplitudes plot
ax = axs[0, 2]
ax.plot(t_axis, amplitudes, color = 'mediumblue')
ax.set_ylabel(r'$A(t)$/V', useTex = True)
ax.set_xlabel(r'$t$/s', useTex = True)
ax.set_title(f"Amplitude", useTex = False)

# spectral of signal plot
ax = axs[1, 0]
f, Pxx_den = periodogram(signal, SAMPLING_F)
print('Spectral Signal')
print(f"Largest/peak frequency found at {f[np.argmax(Pxx_den)]/1e6} MHz")
ax.semilogy(f[1:]/1e3, Pxx_den[1:])
ax.set_ylim(min(Pxx_den[1:])/5, max(Pxx_den[1:])*5)
ax.set_xlabel('frequency [kHz]')
ax.set_ylabel('PSD [V**2/Hz]')
ax.set_title('Periodogram of signal')

# spectral of phase plot
ax = axs[1, 1]
f, Qxx_den = periodogram(phases, SAMPLING_F)
print('Spectral Phase')
print(f"Largest/peak frequency found at {f[np.argmax(Qxx_den)]/1e3} kHz")
ax.semilogy(f[1:]/1e3, Qxx_den[1:])
ax.set_ylim(min(Qxx_den[1:])/5, max(Qxx_den[1:])*5)
ax.set_xlabel('frequency [kHz]')
ax.set_ylabel('PSD [rad**2/Hz]')
ax.set_title('Periodogram of phase')

# spectral of amplitude plot
ax = axs[1, 2]
f, Pxx_den = periodogram(amplitudes, SAMPLING_F)
print('Spectral Amplitude')
print(f"Largest/peak frequency found at {f[np.argmax(Pxx_den)]/1e6} MHz")
ax.semilogy(f[1:]/1e3, Pxx_den[1:])
ax.set_ylim(min(Pxx_den[1:])/5, max(Pxx_den[1:])*5)
ax.set_xlabel('frequency [kHz]')
ax.set_ylabel('PSD [V**2/Hz]')
ax.set_title('Periodogram of amplitudes')



fig.suptitle(f"Mirror $f= ${mirror_f} Hz", useTex= True)
fig.set_size_inches(11.75-1.5, 8.25-2 - 1)
fig.tight_layout()
# fig.subplots_adjust(left=0.10, top=0.91)
plt.show()
plt.pause(15)