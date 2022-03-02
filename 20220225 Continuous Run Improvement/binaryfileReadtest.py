# A test in the ability to read binary files.

# . Honours Module Folder
# ├ FYPLibrary
# | ├ file_reading.py
# | └ IQ_demod.py
# └ Sub Project folder
#   ├ Batch Folders
#   | └ C1-...-000001.txt
#   └ this script(.py)

# Initialiastion: Directory appending for my system. Vary the directories as necessary.
from os import listdir
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
import matplotlib.pyplot as plt
# from matplotlib.ticker import AutoMinorLocator

import lecroyparser as lcp
import binascii

def get_files():
    # uses tkinter to get the paths. returns all files as selected by UI
    # tkinter will not allow for extracting files from multiple folders
    import tkinter as tk
    from tkinter import filedialog

    root = tk.Tk()
    root.withdraw()
    fs = filedialog.askopenfilenames(initialdir= os.path.dirname(__file__), 
        title="Select files")
    return fs

def write_directory(fs):
    # yield folder to write results to.
    # ensure only 1 path. throw error if more than 1
    dir = set(os.path.dirname(f) for f in fs)
    if len(dir) > 1:
        raise ValueError("Files from more than 1 path has been specified.")
    dir = dir.pop()
    # creates write_director (wd) if not already created
    wd = dir + ' Results' # Improper
    my_input = 'Y' # input(f"Results will be written to: {wd}\nAccept? [Y/N]").upper()
    if my_input != 'Y':
        import tkinter as tk
        from tkinter import filedialog
        wd = filedialog.askdirectory()
        print(f"Results will be written to {wd}")
    if not os.path.exists(wd):
        os.makedirs(wd)
    return wd

# Place script specific functions here
from math import log10, floor
def round_to_1(x):
    return round(x, -int(floor(log10(abs(x)))))

# For each file in files selected
def per_file(file, wd, gen_plot, display_plot, **kwargs):
    # file is the full file path
    # wd is the write directory for results of this given file
    # Edit kwargs as necessary
    
    # Read and obtain trace
    NAME = os.path.basename(file)
    if NAME[0] != 'C' and (NAME[-4:] != '.trc' or NAME[-4] != '.txt'):
        print(f"[Warning] Unrecognised file reached, Skipping: {NAME}\n")
        return

    # Edit Initialisation
    SIGNAL_F = 80.625e6*2 #Hz 
    SAMPLING_F = 1.0e7 #Hz
    ph_ad = phase_advance(SIGNAL_F, SAMPLING_F) # phase advance = 2*pi/N
    N, _ = freq_ratio(signal=SIGNAL_F, sample=SAMPLING_F)
    print(f"[Int Log]{N = }")
    # Lecroy parser
    time_axis, signals = parse_and_read_oscilliscope_trc(file)
    signal = signals[0]

    # Get phases over time

    # meta, trace = fr.parse_and_read_oscilliscope_txt(file)
    # [Log] meta = {'Record Length': (500002.0, 'Points'), 'Sample Interval': (1e-06, 's'), 'Trigger Point': (250000.0, 'Samples'), 'Trigger Time': (0.576651, 's'), 'Horizontal Offset': (-0.250001, 's')}
    meta =  {'Record Length': (len(signal), 'Points'), \
        'Sample Interval': (1/SAMPLING_F, 's'), \
        'Trigger Point': ('unknown', 'Samples'), \
        'Trigger Time': ('unknown', 's'), \
        'Horizontal Offset': ('unknown', 's')}
    phases = signal_to_phase(signal, N, ph_ad, phase_advancement_correction= False)
    phases = phase_reconstruction_2(phases, ph_ad)
    t_axis = np.arange(start= 0, 
        stop= (int(meta["Record Length"][0])-N) * meta['Sample Interval'][0], step= meta['Sample Interval'][0])

    # # Do something
    # # <...>
    # value = ...

    # Draw figure
    if gen_plot:
        fig, ax = plt.subplots(nrows=1, ncols=1)
        ax.plot(t_axis, phases, color = 'mediumblue')
        ax.set_ylabel(r'$\phi$/rad', usetex= True)
        ax.set_xlabel(r'$t$/s', usetex= True)
        plt.title(f"Phase Reconstruction", \
            usetex= True )
        if display_plot:
            plt.show(block= False)
            plt.pause(0.8)
        # Saving figure block
        try:
            print('Saving plot...')
            fig.savefig(os.path.join(wd, NAME[:-4] + 'result.png'), 
                format= 'png')
            print('Saved')
        except FileNotFoundError:
            print(f"[Error] Folder should have been created in initialisation block")
            exit()
        plt.close()

    # return value

def final_movement(r):
    # Do something with accumulated values from all traces
    return

def main():
    # print("Select Files to perform this script on.")
    # files = get_files()
    folderpath = os.path.join(os.path.dirname(__file__),'samples')
    files = listdir(path= folderpath)
    print(f"Files selected: {files}")
    DIR_WRITE = write_directory(files)
    print(f"Results will be written to: {DIR_WRITE}")
    results = []
    for file in files:
        results.append(per_file(os.path.join(folderpath, file), 
            DIR_WRITE, True, False))
    final_movement(results)
    print("Script has ended.")

if __name__ == '__main__':
    main()