# This script aims to provide a mathematical model as to how a noise floor
# results in noise in the phase after IQ demodulation

# . Honours Module Folder
# ├ FYPLibrary
# | ├ file_reading.py
# | └ IQ_demod.py
# └ Experiments 1 (Optional Depth)
#   └ Sub Project folder
#     ├ Batch Folders
#     | └ C1-...-000001.txt
#     └ this script(.py)

# Initialiastion: Directory appending for my system. Vary the directories as necessary.
import sys, os.path
# Add the FYP folder for import
if 'FYPLibrary' not in [os.path.basename(x) for x in sys.path]:
    cur_path = os.path.dirname(__file__)
    while 'FYPLibrary' not in os.listdir(cur_path):
        cur_path = os.path.dirname(cur_path)
        if not cur_path:
            raise ValueError('Recursively climbed up directory past root!')
    sys.path.append(os.path.join(cur_path, 'FYPLibrary'))
    del cur_path

# Import Modules
# from os import listdir
from file_reading import *
from IQ_demod import *
import numpy as np
from numpy import random as rd
from numpy import pi
import EPstandard
# from scipy import optimize
from scipy.signal import periodogram
import matplotlib.pyplot as plt
# from matplotlib.ticker import AutoMinorLocator
import matplotlib

# periodogram specific function
# specify sharing axes (https://stackoverflow.com/a/51684195)
def set_share_axes(axs, target=None, sharex=False, sharey=False):
    if target is None:
        target = axs.flat[0]
    # Manage share using grouper objects
    for ax in axs.flat:
        if sharex:
            target._shared_axes['x'].join(target, ax)
        if sharey:
            target._shared_axes['y'].join(target, ax)
    # Turn off x tick labels and offset text for all but the bottom row
    if sharex and axs.ndim > 1:
        for ax in axs[:-1,:].flat:
            ax.xaxis.set_tick_params(which='both', labelbottom=False, labeltop=False)
            ax.xaxis.offsetText.set_visible(False)
    # Turn off y tick labels and offset text for all but the left most column
    if sharey and axs.ndim > 1:
        for ax in axs[:,1:].flat:
            ax.yaxis.set_tick_params(which='both', labelleft=False, labelright=False)
            ax.yaxis.offsetText.set_visible(False)

def main():

    # intialisation
    rd.seed(1337) # seed random
    SIGNAL_F = 80.125e6*2 #Hz 
    SAMPLING_F = 1.0e6 #Hz
    N, _ = freq_ratio(signal=SIGNAL_F, sample=SAMPLING_F)
    ph_ad = phase_advance(SIGNAL_F, SAMPLING_F) # phase advance = 2*pi/N
    print(f"[Int Debug] {N = }, {ph_ad = }")
    amp_V = 0.2 # amplitude of V
    phi_0 = 0.163 # radian associated with t = 0
    mirror_f = 143 # freq of oscillation of mirror
    mirror_A = 0.5 # amplitude of oscillating mirror
    rel_noise = 0.2 # noise of real random relative to ampV
    t_max = 10 #s
    
    t_axis = np.arange(0, t_max, 1/SAMPLING_F)

    signal = amp_V * np.sin(2*pi*SIGNAL_F * t_axis + phi_0 \
        + mirror_A * np.sin(2*pi*mirror_f * t_axis)) \
        + amp_V * rd.default_rng().uniform(-rel_noise, rel_noise, size= len(t_axis))
        # + amp_V * rd.default_rng().normal(0.0, rel_noise, size= len(t_axis))
    phases = signal_to_phase(signal, N, ph_ad, phase_advancement_correction= False)
    phases = phase_reconstruction_2(phases, ph_ad)
    
    # Draw figure
    fig, axs = plt.subplots(nrows= 2, ncols= 2, sharex = False)
    axs[0,0].set_xlim(0, 0.1)
    axs[0,1].set_xlim(0, 0.1)
    set_share_axes(axs[0,:], sharex=True)
    set_share_axes(axs[1, 1:], sharex=True)

    # signal plot
    ax = axs[0, 0]
    try:
        ax.plot(t_axis, signal, color = 'mediumblue')
    except ValueError:
        dt = 1/SAMPLING_F
        axisT = np.fromiter([i* dt for i in range(len(signal))], 
            dtype= np.float64, count= len(signal))
        ax.plot(axisT, signal, color = 'mediumblue')
    ax.set_ylabel(r'$y(t)$/V', usetex= True)
    ax.set_xlabel(r'$t$/s', usetex= True)
    ax.set_title(f"Signal", usetex= False)

    # phase plot
    ax = axs[0, 1]
    ax.plot(t_axis[:1-N], phases, color = 'mediumblue')
    # ax.plot(t_axis, linear(t_axis, *popt), color = 'darkgreen', markersize = 3)
    ax.set_ylabel(r'$\phi$/rad', usetex= True)
    ax.set_xlabel(r'$t$/s', usetex= True)
    ax.set_title(f"Phase", usetex= False)

    # spectral of signal plot
    ax = axs[1, 0]
    f, power_den = periodogram(signal, SAMPLING_F)
    print('\nSpectral Signal')
    print(f"Largest/peak frequency found at {f[np.argmax(power_den)]/1e3} kHz")
    ax.semilogy(f[0:]/1e3, power_den[0:])
    # ax.set_xscale('log')
    ax.set_ylim(min(power_den[1:])/5, max(power_den[1:])*5)
    ax.set_xlabel('frequency [kHz]')
    ax.set_ylabel('PSD [V**2/Hz]')
    ax.set_title('Signal FFT')

    # spectral of phase plot
    ax = axs[1, 1]
    f, Qxx_den = periodogram(phases, SAMPLING_F)
    print('\nSpectral Phase')
    print(f"Largest/peak frequency found at {f[np.argmax(Qxx_den)]/1e3} kHz")
    ax.semilogy(f[0:]/1e3, Qxx_den[0:])
    # ax.set_xscale('log')
    ax.set_ylim(min(Qxx_den[1:])/5, max(Qxx_den[1:])*5)
    ax.set_xlabel('frequency [kHz]')
    ax.set_ylabel('PSD [rad**2/Hz]')
    ax.set_title(f'phase FFT\nPeak at {f[np.argmax(Qxx_den)]/1e3 :.6f} kHz')

    fig.set_size_inches(11.75-1.5, 8.25-2 - 1)
    fig.tight_layout()

    if True:
        plt.show(block= True)
        plt.pause(0.8)

    return


if __name__ == '__main__':
    main()
