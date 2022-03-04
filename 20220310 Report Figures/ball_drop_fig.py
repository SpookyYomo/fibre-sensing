# Generate ball drop from 20210914 Initial Ball Drop\C2-ball-drop-00002.txt

# . Honours Module Folder
# ├ FYPLibrary
# | ├ file_reading.py
# | └ IQ_demod.py
# ├ 20210914 Initial Ball Drop
# | └ C2-ball-drop-00002.txt
# └ Sub Project folder
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
    my_input = input(f"Results will be written to: {wd}\nAccept? [Y/N]").upper()
    if my_input != 'Y':
        import tkinter as tk
        from tkinter import filedialog
        wd = filedialog.askdirectory()
        print(f"Results will be written to {wd}")
    if not os.path.exists(wd):
        os.makedirs(wd)
    return wd

def main():
    FYPLib = os.path.dirname(os.path.dirname(__file__))
    file = os.path.join(FYPLib, 
        "20210914 Initial Ball Drop", "C2-ball-drop-00002.txt")
    print(f"File selected: {file}")
    wd = write_directory(__file__)
    print(f"Results will be written to: {wd}")

    # Edit Initialisation
    SIGNAL_F = 80.0125e6*2 #Hz 
    SAMPLING_F = 1.0e5 #Hz
    ph_ad = phase_advance(SIGNAL_F, SAMPLING_F) # phase advance = 2*pi/N
    N, _ = freq_ratio(signal=SIGNAL_F, sample=SAMPLING_F)
    print(f"[Int Debug] {N = }, {ph_ad*2/np.pi = } should be 1")

    readTxt = read_oscilliscope_txt(file)
    meta, trace = parse_oscilliscope_txt(readTxt)
    signal = signal_from_trace(trace)
    phases = signal_to_phase(signal, N, ph_ad, phase_advancement_correction= False)
    phases = phase_reconstruction_2(phases, ph_ad)
    t_axis = np.arange(start= 0, 
        stop= (int(meta["Record Length"][0])-N+1) * meta['Sample Interval'][0], step= meta['Sample Interval'][0])

    matplotlib.rcParams.update({'font.size': 11})

    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.plot(t_axis, phases, color = 'mediumblue')
    ax.set_ylabel(r'$\phi_d$/rad', usetex= True)
    ax.set_xlabel(r'$t$/s', usetex= True)
    ax.set_ylim([-50, 40])
    ax.set_xlim([3.2, 4.0])

    plt.title(f"Table tap, $f_s$ = {SAMPLING_F/1e3:.1f} kS/s, \
        $2\Omega$ = {SIGNAL_F/1e6:.3f} MHz, $N$ = {N}", usetex= True )

    # powerpoint is 13.333 inches wide by 7.5 inches high
    # geometry package boundary -0.875in*2 horizontally
    # figures set to 0.85 * textwidth
    fig.set_size_inches(0.85*(8.25-0.875*2), 3.6) 
    fig.tight_layout()

    # Saving figure block
    try:
        print('Saving plot...')
        NAME = os.path.basename(file)
        fig.savefig(os.path.join(wd, NAME[:-4] + 'result.png'), 
            format= 'png', dpi= 300, pad_inches = 0.01)
        print('Saved')
    except FileNotFoundError:
        print(f"[Error] Folder should have been created in initialisation block")
        exit()
    plt.close()

if __name__ == '__main__':
    main()