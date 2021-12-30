# This script aims to help identify the phi(t) with different sampling frequencies. 
# You will need this script to be placed along side the folder 
# containing the various .txt files of the oscilliscopes that you wish 
# to identify. 

# . Honours Module Folder
# ├ FYPLibrary
# | ├ file_reading.py
# | └ IQ_demod.py
# └ Sub Project folder
#   ├ this script(.py)
#   └ "Batch" folder
#     └ C2-xxxHz-xxxmVpp-80.125xxxMHz-comment.txt

# Oscilliscopes reading are to contain header data, time series and 
# values, and are space seperated

# Initialiastion: Directory appending for my system. Vary the directories as necessary.
import sys, os.path
# Add the FYPLibrary folder in the same level as this folder.
if os.path.dirname(os.path.dirname(__file__)) not in sys.path:
    # sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'FYPLibrary'))

# from pprint import pprint
# pprint(sys.path)

# Import Modules
from os import listdir
from file_reading import *
from IQ_demod import *
import numpy as np
import matplotlib.pyplot as plt
# from EPstandard import easy_read_popt_pcov
# from scipy import optimize
# from matplotlib.ticker import AutoMinorLocator

# Step 1: Change the name to which the folder you wish to read from
# Yield all files in batch folder
folderpath = os.path.join(os.path.dirname(__file__),'batch 1')

def main():
    # Step 2: Change functionalities you wish to perform
    INDIVIDUAL_PLOTS = True
    DISPLAY_INDIVIDUAL_PLOTS = False

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

    # Plot for every trace
    for file in files:
        NAME = file
        TXT_FILE_PATH = file_str_to_path(NAME, dir= folderpath)
        if NAME[0] != 'C' and NAME[-4:] != '.txt':
            print(f"[Warning] Unrecognised file reached, Skipping: {NAME}\n")
            continue

        _, f, mVpp, SAMPLING_F = numerical_variables_from_name(NAME)
        SAMPLING_F *= 1e6
        SIGNAL_F = 160e6 + SAMPLING_F / 4
        N, _ = freq_ratio(signal=SIGNAL_F, sample=SAMPLING_F)
        print(f"[Int Debug] {NAME = }, {f = }, {mVpp = }, {SAMPLING_F = }, {N = }")

        meta, trace = fr.parse_and_read_oscilliscope_txt(TXT_FILE_PATH)
        signal = signal_from_trace(np.asarray(trace))
        phases = signal_to_phase(signal, N, 2*pi/N, phase_advancement_correction=True)
        phases = phase_reconstruction(phases, 4.3)
        # phases = signal_to_phase(signal, N, 2*pi/N, phase_advancement_correction=False)
        # phases = phase_reconstruction_2(phases, 2*pi/N)
        t_axis = np.arange(start= 0, 
            stop= (int(meta["Record Length"][0])-N) * meta['Sample Interval'][0], step= meta['Sample Interval'][0])

        if INDIVIDUAL_PLOTS:
            fig, ax = plt.subplots(nrows=1, ncols=1)
            ax.plot(t_axis, phases, color = 'mediumblue')
            ax.set_ylabel(r'$\phi$/rad', useTex = True)
            ax.set_xlabel(r'$t$/s', useTex = True)
            plt.title(f"Phase Reconstruction 1:$f_{{Sampling}}$ = {SAMPLING_F/1e6}MHz\n$f_{{Mirror}}$ = {f}Hz, Piezo Ampl = {mVpp}mVpp", \
                useTex= True )
            if DISPLAY_INDIVIDUAL_PLOTS:
                plt.show(block= False)
                plt.pause(0.8)
            # Saving figure block
            try:
                print('Saving plot...')
                fig.savefig(os.path.join(DIRECTORY_WRITE, NAME[:-4] + 'method1.png'), 
                    format= 'png')
                print('Saved')
            except FileNotFoundError:
                print(f"[Error] Folder should have been created in initialisation block")
                exit()
            plt.close()

if __name__ == '__main__':
    main()