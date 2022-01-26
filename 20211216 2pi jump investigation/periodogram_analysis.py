# This script aims to get the power spectrum for multiple sampling
# frequencies, those without and without table tapping occuring. Due to 
# prior investigation, it has been noted that phase_reconstruction_2 is 
# preferred, particularly for sampling rate of 1MHz.

# . Honours Module Folder
# ├ FYPLibrary
# | ├ file_reading.py
# | └ IQ_demod.py
# └ Sub Project folder
#   ├ Batch Folders
#   | └ C1-...-000001.txt
#   └ this script(.py)

# Initialiastion: Directory appending for my system. Vary the directories as necessary.
import sys, os.path

# from numpy.core.fromnumeric import argmax
# Add the FYPLibrary folder in the same level as this folder.
if os.path.dirname(os.path.dirname(__file__)) not in sys.path:
    # print(f"{os.path.dirname(os.path.dirname(__file__)) = }")
    # sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'FYPLibrary'))

# from pprint import pprint
# pprint(sys.path)

# Import Modules
from os import listdir
from file_reading import *
from IQ_demod import *
import numpy as np
from numpy import float32, pi as pi
import EPstandard
from scipy.signal import periodogram
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# periodogram specific function
# specify sharing axes
def set_share_axes(axs, target=None, sharex=False, sharey=False):
    if target is None:
        target = axs.flat[0]
    # Manage share using grouper objects
    for ax in axs.flat:
        if sharex:
            target._shared_x_axes.join(target, ax)
        if sharey:
            target._shared_y_axes.join(target, ax)
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
    # Step 1: Change the according folder
    # Yield file in batch folder
    folderpath = os.path.dirname(__file__)+'\Batch 1'

    # Step 2: Change functionalities you wish to perform
    DISPLAY_INDIVIDUAL_PLOTS = True

    # Initialisation
    # results = [] # (frequency deviation, gradient (m), std deviation (delta_m)) 
    DIRECTORY_WRITE = folderpath + ' Results' # Improper
    if not os.path.exists(DIRECTORY_WRITE):
        os.makedirs(DIRECTORY_WRITE)
    files = listdir(path= folderpath)

    for i, f in enumerate(files):
        print(f"{i:>3}: {f}")
    print("Select -1 for all")
    my_input = int(input("Please specify which file you would like to generate a periodogram of: "))
    if my_input != -1:
        files = (files[my_input], )

    for file in files:
        NAME = file
        TXT_FILE_PATH = file_str_to_path(NAME, dir=folderpath)
        if NAME[0] != 'C' and NAME[-4:] != '.txt':
            print(f"[Warning] Unrecognised file reached, Skipping: {NAME}\n")
            continue

        _, mirror_f, mVpp, SAMPLING_F = numerical_variables_from_name(NAME)
        SAMPLING_F *= 1e6
        SIGNAL_F = 160e6 + SAMPLING_F / 4
        N, _ = freq_ratio(signal=SIGNAL_F, sample=SAMPLING_F)
        print(f"[Int Debug] {NAME = }, {mirror_f = }, {mVpp = }, {SAMPLING_F = }, {N = }")

        # file reading
        trace_meta, trace = fr.parse_and_read_oscilliscope_txt(TXT_FILE_PATH)
        signal = signal_from_trace(np.asarray(trace))
        # phases = signal_to_phase(signal, N, 2*pi/N, phase_advancement_correction=True)
        # phases = phase_reconstruction(phases, 4.3)
        phases = signal_to_phase(signal, N, 2*pi/N, phase_advancement_correction=False)
        phases = phase_reconstruction_2(phases, 2*pi/N)
        amplitudes = get_R_signal(signal, N, 2*pi/N) # amplitude over time
        t_axis = np.arange(start= 0, 
            stop= (int(trace_meta["Record Length"][0])-N) * trace_meta['Sample Interval'][0], 
            step= trace_meta['Sample Interval'][0])
        time_axis = np.arange(start= 0, 
            stop= trace_meta["Record Length"][0] * trace_meta['Sample Interval'][0], 
            step= trace_meta['Sample Interval'][0])

        # generate plot
        fig, axs = plt.subplots(nrows = 2, ncols= 3, sharex= False)
        
        set_share_axes(axs[0,:], sharex=True)
        set_share_axes(axs[1,:], sharex=True)

        # signal plot
        ax = axs[0, 0]
        try:
            ax.plot(time_axis, signal, color = 'mediumblue')
        except ValueError:
            dt = trace_meta['Sample Interval'][0]
            time_axis = np.fromiter([i* dt for i in range(len(signal))], 
                dtype= np.float32, count= len(signal))
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

        print()
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
        ax.set_title(f'Periodogram of phase\nPeak at {f[np.argmax(Qxx_den)]/1e3 :.6f} kHz')

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


        fig.suptitle(f"Sampling freq $f= ${SAMPLING_F/1e6} MHz, AOM freq = f{SIGNAL_F/1e6:.6f}MHz", useTex= True)
        fig.set_size_inches(11.75-1.5, 8.25-2 - 1)
        fig.tight_layout()
        # fig.subplots_adjust(left=0.10, top=0.91)
        if DISPLAY_INDIVIDUAL_PLOTS:
            plt.show(block = (my_input!=-1))
            plt.pause(15)
        print("Saving fig...")
        fig.savefig(os.path.join(DIRECTORY_WRITE, NAME[:-4] + '_periodogram.png'), \
            format= 'png')
        print("Saved\n")
    print('All plots have been completed')

if __name__ == '__main__':
    main()