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
# Add the FYP folder for import
if 'FYPLibrary' not in [os.path.basename(x) for x in sys.path]:
    cur_path = os.path.dirname(__file__)
    while 'FYPLibrary' not in os.listdir(cur_path):
        cur_path = os.path.dirname(cur_path)
        if not cur_path:
            raise ValueError('Recursively climbed up directory past root!')
    sys.path.append(os.path.join(cur_path, 'FYPLibrary'))
    del cur_path

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
    import sys
    from PyQt5 import QtWidgets

    class get_files(QtWidgets.QMainWindow):
        def __init__(self):
            super().__init__()
            self._filename, _ = QtWidgets.QFileDialog.getOpenFileNames()
        
        @property
        def filename(self):
            return self._filename

    # generates app first for window to populate
    app = QtWidgets.QApplication(sys.argv)
    fs = get_files().filename        
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
    my_input = input(f"Results will be written to: {wd}\nAccept? [Y/N]")
    if my_input == 'N':
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
    if NAME[0] != 'C' and NAME[-4:] != '.txt':
        print(f"[Warning] Unrecognised file reached, Skipping: {NAME}\n")
        return

    # Edit Initialisation
    num_vars = numerical_variables_from_name(NAME)
    SIGNAL_F = 80.125e6*2 #Hz 
    SAMPLING_F = 1.0e6 #Hz
    ph_ad = phase_advance(SIGNAL_F, SAMPLING_F) # phase advance = 2*pi/N
    N, _ = freq_ratio(signal=SIGNAL_F, sample=SAMPLING_F)
    print(f"[Int Debug] {num_vars = } {N = }")
    
    # Get phases over time
    if file[-4:] == '.txt':
        meta, trace = parse_and_read_oscilliscope_txt(file)
        signal = signal_from_trace(np.asarray(trace))
    elif file[-4:] == '.trc':
        time, [signal,] = parse_and_read_oscilliscope_trc(file)
        meta =  {'Record Length': (len(signal), 'Points'), \
            'Sample Interval': (1/SAMPLING_F, 's'), \
            'Trigger Point': ('unknown', 'Samples'), \
            'Trigger Time': ('unknown', 's'), \
            'Horizontal Offset': ('unknown', 's')}
    else:
        raise ValueError("Wrong file type passed.")
    phases = signal_to_phase(signal, N, ph_ad, phase_advancement_correction= False)
    phases = phase_reconstruction_2(phases, ph_ad)
    t_axis = np.arange(start= 0, 
        stop= (int(meta["Record Length"][0])-N+1) * meta['Sample Interval'][0], step= meta['Sample Interval'][0])

    # Do something
    # <...>

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

    return

def main():
    print("Select Files to perform this script on.")
    files = get_files()
    print(f"Files selected: {files}")
    DIR_WRITE = write_directory(files)
    print(f"Results will be written to: {DIR_WRITE}")
    for file in files:
        per_file(file, DIR_WRITE, True, False)
    print("Script has ended.")

if __name__ == '__main__':
    main()