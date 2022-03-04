# Intended to work on binary trc files that is reading continuous data 
# over large periods of time. Writes to log: 
# (id, mean, std_dev, max, min)

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
# import matplotlib.pyplot as plt
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
    # yield folder to write results to.
    # ensure only 1 path. throw error if more than 1
    dir = set(os.path.dirname(f) for f in fs)
    if len(dir) > 1:
        raise ValueError("Files from more than 1 path has been specified.")
    dir = dir.pop()
    # creates write_director (wd) if not already created
    wd = dir + ' Results' # Improper
    my_input = input(f"Results will be written to: {wd}\nAccept? [Y/N] ")
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
def per_file(file, logFile):
    # log file for writing 
    
    # Read and obtain trace
    NAME = os.path.basename(file)
    if NAME[0] != 'C' and NAME[-4:] != '.trc':
        print(f"[Warning] Unrecognised file reached, Skipping: {NAME}\n")
        return

    # Edit Initialisation
    id = NAME[-9:-4]
    SIGNAL_F = 70.625e6*2 #Hz 
    SAMPLING_F = 1.0e7 #Hz
    ph_ad = phase_advance(SIGNAL_F, SAMPLING_F) # phase advance = 2*pi/N
    N, _ = freq_ratio(signal=SIGNAL_F, sample=SAMPLING_F)
    print(f"[Int Debug] {NAME = }", f"{N = }", f"{id = }", sep= "[Int Debug]")
    
    # Lecroy read
    _, signal = fr.parse_and_read_oscilliscope_trc(file)
    signal = signal[0]
    phases = signal_to_phase(signal, N, ph_ad, phase_advancement_correction= False)
    phases = phase_reconstruction_2(phases, ph_ad)

    logLine = [id, np.mean(phases), np.std(phases), np.max(phases), np.min(phases)]
    logLine = np.asarray(logLine, dtype= str)

    with open(logFile, "a+") as file_object:
        # Move read cursor to the start of file.
        file_object.seek(0)
        # If file is not empty then append '\n'
        if len(file_object.read(100)) > 0 :
            file_object.write("\n")
        # Append text at the end of file
        file_object.write("\t".join(logLine))

    return

def main():
    print("Select Files to perform this script on.")
    files = get_files()
    for f in files:
        NAME = os.path.basename(f)
        if NAME[0] != 'C' and NAME[-4:] != '.trc':
            print(f"[Warning] Unrecognised file reached, Skipping: {NAME}\n")
            return
    print(f"Number of files selected: {len(files)}")

    print("Select log to save processed data to.")
    log_file = get_files()
    log_file = log_file[0]
    if log_file[-4:] != '.txt':
        print(f"[Warning] Selected file being {log_file} cannot be logged to!")
        return
    print(f"Results will be written to: {log_file}")
    for file in files:
        per_file(file, log_file)
    print("Script has ended.")

if __name__ == '__main__':
    main()