# Base Script to edit off from. Acts on a trace file.

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

# Add the FYPLibrary folder in the same level as this folder.
if os.path.dirname(os.path.dirname(__file__)) not in sys.path:
    sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'FYPLibrary'))

# from pprint import pprint
# pprint(sys.path)

# Import Modules
# from os import listdir
from file_reading import *
from IQ_demod import *
import numpy as np
# from numpy import pi
# import EPstandard
# from scipy.signal import periodogram
import matplotlib.pyplot as plt
# from matplotlib.ticker import AutoMinorLocator

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

# Place script specific functions here
def my_func():
    return

# For each file in files selected
def per_file(file, wd, gen_plot, display_plot, **kwargs):
    # wd is the write directory for results of this given file
    # Edit kwargs as necessary
    
    # Read and obtain trace
    NAME = os.path.basename(file)
    if NAME[0] != 'C' and NAME[-4:] != '.trc':
        print(f"[Warning] Unrecognised file reached, Skipping: {NAME}\n")
        return

    # Edit Initialisation
    _, SAMPLING_F = numerical_variables_from_name(NAME)
    SAMPLING_F *= 1e6
    SIGNAL_F = 140e6 + SAMPLING_F/4 #Hz 
    ph_ad = phase_advance(SIGNAL_F, SAMPLING_F) # phase advance = 2*pi/N
    N, _ = freq_ratio(signal=SIGNAL_F, sample=SAMPLING_F)
    print(f"[Int Debug] {N = }, {ph_ad*2/np.pi = } should be 1")
    print(f"[Int Debug] {NAME = }, SAMPLING_F = {SAMPLING_F/1e6} MHz")
    
    # Lecroy parser
    time_axis, signals = parse_and_read_oscilliscope_trc(file)
    signal = signals[0]
    meta =  {'Record Length': (len(signal), 'Points'), \
        'Sample Interval': (1/SAMPLING_F, 's'), \
        'Trigger Point': ('unknown', 'Samples'), \
        'Trigger Time': ('unknown', 's'), \
        'Horizontal Offset': ('unknown', 's')}
    phases = signal_to_phase(signal, N, ph_ad, phase_advancement_correction= False)
    deltas = phase_reconstruction_2(phases, ph_ad, get_deltas=True)
    phases = phase_reconstruction_2(phases, ph_ad)
    t_axis = np.arange(start= 0, 
        stop= (int(meta["Record Length"][0])-N+1) * meta['Sample Interval'][0], step= meta['Sample Interval'][0])

    # Do something
    value = [str(SAMPLING_F/1e6), str(np.std(deltas))[:7]]

    # Draw figure
    if gen_plot and SAMPLING_F == 1e6:
        fig, ax = plt.subplots(nrows=1, ncols=1)
        ax.plot(t_axis, phases, color = 'mediumblue', linewidth= 0.1)
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
            fig.savefig(os.path.join(wd, NAME[:-4] + '_phi_d.png'), 
                format= 'png', dpi= 300, pad_inches= 0)
            print('Saved')
        except FileNotFoundError:
            print(f"[Error] Folder should have been created in initialisation block")
            exit()
        plt.close()

        fig, ax = plt.subplots(nrows=1, ncols=1)
        ax.hist(deltas, bins= 615)
        ax.set_yscale('log')
        ax.set_xlim(np.min(deltas)-0.1, np.max(deltas)+0.1)
        ax.set_ylabel(r'Counts', usetex= False)
        ax.set_xlabel(r'$\Delta\phi$/rad', usetex= True)
        # figures set to 0.59 * textwidth # DIFFERENT!!
        fig.set_size_inches(0.59*(8.25-1.0*2), 1.8) 
        fig.tight_layout()
        if display_plot:
            plt.show(block= False)
            plt.pause(0.8)
        # Saving figure block
        try:
            print('Saving plot...')
            fig.savefig(os.path.join(wd, 'hist_' + NAME[:-4] + '_hist.png'), 
                format= 'png', dpi= 300, pad_inches= 0)
            print('Saved')
        except FileNotFoundError:
            print(f"[Error] Folder should have been created in initialisation block")
            exit()
        plt.close()

    return value

def final_movement(rv):
    # Do something with accumulated values from all traces
    for i in rv:
        print(" & ".join(i), end= " \\\\ \n")
    print(r"\hline")
    return

def main():
    print("Select Files to perform this script on.")
    files_folder = os.path.join(os.path.dirname(os.path.dirname(__file__)), 
        "20220316 130 meter cable", "Batch 1")
    files = [os.path.join(files_folder,f) for f in os.listdir(files_folder) if os.path.isfile(os.path.join(files_folder,f))]
    print(f"Files selected: {files}")
    DIR_WRITE = os.path.join(write_directory(__file__), "phase_distribution")
    print(f"Results will be written to: {DIR_WRITE}")
    results = []
    for file in files:
        results.append(per_file(file, DIR_WRITE, True, False))
    final_movement(results)
    print("Script has ended.")

if __name__ == '__main__':
    main()