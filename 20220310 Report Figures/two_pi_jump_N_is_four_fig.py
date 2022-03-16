# Generate 2 pi jump when N = 4 from 20211214 Linearity with Amp and Freq Incomplete\Batch1\C1-004Hz-020.4mVpp-80.125000MHz-OneAmp00001.txt

# . Honours Module Folder
# ├ FYPLibrary
# | ├ file_reading.py
# | └ IQ_demod.py
# └ Sub Project folder
#   ├ traces folder
#   | └ C1-004Hz-020.4mVpp-80.125000MHz-OneAmp00001.txt
#   └ this script(.py)

# Initialiastion: Directory appending for my system. Vary the directories as necessary.
import sys, os.path

# Add the FYPLibrary folder in the same level as this folder.
if os.path.dirname(os.path.dirname(__file__)) not in sys.path:
    sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'FYPLibrary'))

# from pprint import pprint
# pprint(sys.path)

# Import Modules
from os import listdir
from file_reading import *
from IQ_demod import *
import numpy as np
# from numpy import pi
# import EPstandard
# from scipy.signal import periodogram
import matplotlib
import matplotlib.pyplot as plt
# from matplotlib.ticker import AutoMinorLocator
import string

def write_directory(fs):
    # modified from the usual write dir
    # yield folder to write results to.
    # ensure only 1 path. throw error if more than 1
    dir = os.path.dirname(fs)
    # if len(dir) > 1:
    #     raise ValueError("Files from more than 1 path has been specified.")
    # dir = dir.pop()
    # creates write_director (wd) if not already created
    wd = dir # Improper
    # my_input = input(f"Results will be written to: {wd}\nAccept? [Y/N]").upper()
    my_input = 'Y'
    if my_input != 'Y':
        import tkinter as tk
        from tkinter import filedialog
        wd = filedialog.askdirectory()
        print(f"Results will be written to {wd}")
    if not os.path.exists(wd):
        os.makedirs(wd)
    return wd

def main():
    file = os.path.join(os.path.dirname(__file__), 
        "traces", "C1-004Hz-020.4mVpp-80.125000MHz-OneAmp00001.txt")
    print(f"File selected: {file}")
    wd = write_directory(__file__)
    print(f"Results will be written to: {wd}")

    # Edit Initialisation
    SIGNAL_F = 80.125e6*2 #Hz 
    SAMPLING_F = 1.0e6 #Hz
    ph_ad = phase_advance(SIGNAL_F, SAMPLING_F) # phase advance = 2*pi/N
    N, _ = freq_ratio(signal=SIGNAL_F, sample=SAMPLING_F)
    print(f"[Int Debug] {N = }, {ph_ad*2/np.pi = } should be 1")

    meta, trace = parse_and_read_oscilliscope_txt(file)
    signal = signal_from_trace(trace)
    phases = signal_to_phase(signal, N, ph_ad, phase_advancement_correction= False)
    phases = phase_reconstruction_2(phases, ph_ad)
    t_axis = np.arange(start= 0, 
        stop= (int(meta["Record Length"][0])-N+1) * meta['Sample Interval'][0], step= meta['Sample Interval'][0])
    radii = get_R_signal(signal, N, ph_ad)

    matplotlib.rcParams.update({'font.size': 11})

    fig, axs = plt.subplots(nrows=2, ncols=1, sharex= True)
    ax = axs[0]
    ax.plot(t_axis, phases, color = 'mediumblue')
    ax.set_ylabel(r'$\phi_d$/rad', usetex= True)
    # ax.set_xlabel(r'$t$/s', usetex= True)
    ax.set_ylim([min(phases) - 0.5, max(phases) + 0.5])
    ax.set_xlim([0.0, 0.5])
    ax.axes.xaxis.set_visible(False)

    ax = axs[1]
    ax.plot(t_axis, radii, color = 'mediumblue')
    ax.set_ylabel(r'$r$', usetex= True)
    ax.set_xlabel(r'$t$/s', usetex= True)
    ax.set_xlim([0.0, 0.5])

    for n, ax in enumerate(axs.flat):  
        ax.text(-0.1, 1.1, "(" + string.ascii_lowercase[n] + ")", transform=ax.transAxes, 
                size=11, weight='bold')

    # powerpoint is 13.333 inches wide by 7.5 inches high
    # geometry package boundary -0.875in*2 horizontally
    # figures set to 0.85 * textwidth
    fig.set_size_inches(0.85*(8.25-0.875*2), 3.6) 
    fig.tight_layout()

    # Saving figure block
    try:
        print('Saving plot...')
        NAME = os.path.basename(file)
        fig.savefig(os.path.join(wd, NAME[:-4] + '_2pijump.png'), 
            format= 'png', dpi= 300, pad_inches = 0.01)
        print('Saved')
    except FileNotFoundError:
        print(f"[Error] Folder should have been created in initialisation block")
        exit()
    plt.close()

if __name__ == '__main__':
    main()