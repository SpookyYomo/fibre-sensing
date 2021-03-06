# This script aims to debug the NEW (array-level) operation of 
# phase_reconstruction_2, i.e.: phase accumulation model. The choice of 
# trace is given by an explicit path.

# . Honours Module Folder
# ├ FYPLibrary
# | ├ file_reading.py
# | └ IQ_demod.py
# └ Experiments 1 (Optional Depth)
#   └ Sub Project folder
#     └ this script(.py)

# Initialization: Directory appending for my system. Vary the directories as necessary.
import sys, os
# Add the FYPLibrary folder in the same level as this folder.
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
from os import listdir
from file_reading import *
from IQ_demod import *
import numpy as np
import matplotlib.pyplot as plt
from math import pi
# from EPstandard import easy_read_popt_pcov
# from scipy import optimize
# from matplotlib.ticker import AutoMinorLocator

def main():
    folderpath = r'D:\NUS\Semester 7\PC4199 Honours Project\20211215 phase reconstruction experimental\traces'
    for i, f in enumerate(files:= listdir(path= folderpath)):
        print(f"{i:>3}: {f}")
    # my_input = int(input("Please specify which file you would like to generate a periodogram of: "))
    my_input = 5
    NAME = files[my_input]
    TXT_FILE_PATH = file_str_to_path(NAME, dir=folderpath)
    print(f"{TXT_FILE_PATH = }")

    SIGNAL_F = 80.125e6 *2
    SAMPLING_F = 1e6
    N, _ = freq_ratio(signal= SIGNAL_F, sample= SAMPLING_F)

    print(f"{N = }")

    meta, trace = parse_and_read_oscilliscope_txt(TXT_FILE_PATH)
    signal = signal_from_trace(np.asarray(trace))
    phases = signal_to_phase(signal, N, 2*pi/N, phase_advancement_correction= False)
    
    phases = phase_reconstruction_2(phases, phase_advance= 2*pi/N)
    s=8
    phases_avg = np.convolve(phases, np.ones(s), "valid")/s
    t_axis = np.arange(start= 0, 
        stop= (int(meta["Record Length"][0])-N+1) * meta['Sample Interval'][0],
        step= meta['Sample Interval'][0])[:-s+1]

    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.plot(t_axis, phases_avg, color = 'mediumblue')
    ax.set_ylabel(r'$\phi$/rad', useTex = True)
    ax.set_xlabel(r'$t$/s', useTex = True)
    fig.set_size_inches(11.75-1.5, 8.25-2 - 1)
    plt.show(block= True)

if __name__ == '__main__':
    main()