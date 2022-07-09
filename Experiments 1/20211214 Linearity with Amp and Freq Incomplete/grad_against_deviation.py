# This script aims to help identify the correction to the laser frequency 
# that allows for the least drift due to a mismatch in the driving frequency 
# and the oscilliscope's internal clock. You will need this script to be 
# placed along side the folder containing the various .txt files of the 
# oscilliscopes that you wish to identify. This script is for *non-zero*
# mirror frequencies.

# . Honours Module Folder
# ├ FYPLibrary
# | ├ file_reading.py
# | └ IQ_demod.py
# └ Experiments 1 (Optional Depth)
#   └ Sub Project folder
#     ├ this script(.py)
#     └ "Batch" folder
#       └ C2-xxxHz-xxxmVpp-80.125xxxMHz-comment.txt

# Oscilliscopes reading are to contain header data, time series and 
# values, and are space seperated

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
from os import listdir
from file_reading import *
from IQ_demod import *
import numpy as np
import matplotlib.pyplot as plt
import EPstandard
from scipy import optimize
from matplotlib.ticker import AutoMinorLocator

# Step 1: Change the name to which the folder you wish to read from
# Yield all files in batch folder
folderpath = os.path.dirname(__file__)+'\Batch 0'


# This file's specific function
def linear_root(x, m, r):
    return m*(x-r)

def drifting_sin(x, A, f, phi, m, c):
    return A * np.sin(2*pi*f*x + phi) + m*x+c

def main():
    # Step 2: Change functionalities you wish to perform
    INDIVIDUAL_PLOTS = True
    DISPLAY_INDIVIDUAL_PLOTS = False
    FINAL_REGRESSION = True

    # Initialisation
    SIGNAL_F = 80.125e6*2 #Hz 
    SAMPLING_F = 1.0e6 #Hz
    N, _ = freq_ratio(signal=SIGNAL_F, sample=SAMPLING_F)
    results = [] # (frequency deviation, gradient (m), std deviation (delta_m)) 
    DIRECTORY_WRITE = folderpath + ' Results' # Improper
    if not os.path.exists(DIRECTORY_WRITE):
        os.makedirs(DIRECTORY_WRITE)
    files = listdir(path= folderpath)
    
    # Section A: Determining the drifting for each trace.
    for file in files:
        NAME = file
        TXT_FILE_PATH = file_str_to_path(NAME, dir=folderpath)
        if NAME[0] != 'C':
            print(f"[Warning] Unrecognised file reached, Skipping: {NAME}\n")
            continue
        
        # delta corr should have the unit of Hz here.
        _, f, mVpp, delta_corr = numerical_variables_from_name(NAME)
        if f == 0:
            print(F"[Warning] This script requires a mirror driving frequency! The file will be skipped: {NAME}")
            continue
        delta_corr -= 80.125
        delta_corr *= 1e6 # converts MHz to Hz
        print(f"[Int Debug] {NAME = }, {f = }, {mVpp = }, {delta_corr = }")
        
        meta, trace = fr.parse_and_read_oscilliscope_txt(TXT_FILE_PATH)
        signal = signal_from_trace(np.asarray(trace))
        phases = signal_to_phase(signal, N, 2*pi/N, True)
        phases = phase_reconstruction(phases, 4.3)
        t_axis = np.arange(start= 0, 
            stop= (int(meta["Record Length"][0])-N+1) * meta['Sample Interval'][0], step= meta['Sample Interval'][0])
        
        drifting_sin_bounds = ((0.03*mVpp,      0,  -2*pi, -np.inf, -np.inf), \
                               (   np.inf, np.inf,  2*pi,  np.inf,  np.inf))
        popt, pcov = optimize.curve_fit(drifting_sin, t_axis, phases, 
            p0= [0.05*mVpp, f, 0, delta_corr*12.4, 0.5*(max(phases[:1000])+min(phases[:1000]))], bounds= drifting_sin_bounds)
        fittings = EPstandard.easy_read_popt_pcov(popt, pcov)
        
        if INDIVIDUAL_PLOTS:
            fig, ax = plt.subplots(nrows=1, ncols=1)
            ax.plot(t_axis, phases, color = 'mediumblue')
            ax.plot(t_axis, drifting_sin(t_axis, *popt), color = 'darkgreen', markersize = 3)
            ax.set_ylabel(r'$\phi$/rad', useTex = True)
            ax.set_xlabel(r'$t$/s', useTex = True)
            plt.title(f"$f_{{Mirror}}$ = {f}Hz, Piezo Ampl = {mVpp}mVpp, \
                $\Delta$ = {delta_corr:.5f}Hz, \n(A, $f$, $\phi_0$, $m$, \
                    c) = {fittings}", useTex= True )
            if DISPLAY_INDIVIDUAL_PLOTS:
                plt.show(block= False)
                plt.pause(0.8)
            # Saving figure block
            try:
                print('Saving plot...')
                fig.savefig(os.path.join(DIRECTORY_WRITE, NAME[:-4] + '.png'), 
                    format= 'png')
                print('Saved')
            except FileNotFoundError:
                print(f"[Error] Folder should have been created in initialisation block")
                exit()
            plt.close()
        
        results.append((delta_corr, popt[3], np.sqrt(pcov[3][3])))
        print()

    xs = [r[0] for r in results]
    ys = [r[1] for r in results]
    yerrs = [r[2] for r in results]

    del fittings # prepare for final regression component

    # Section B: Plot of drift against laser correction frequencies.
    # final plot of all gradients against various frequencies
    print("Commencing final diagram...")
    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.errorbar(xs, ys, 
        yerr=yerrs, fmt = 'bo', markersize = 3)

    if FINAL_REGRESSION:
        popt, pcov = optimize.curve_fit(linear_root, xs, ys, p0=(100, 55), sigma= yerrs)
        fittings = EPstandard.easy_read_popt_pcov(popt,pcov)
        ax.plot(delta_axis:= np.arange(start=np.min(xs), stop=np.max(xs), step = 1), 
            linear_root(delta_axis, *popt),
            markersize = 1)
        plt.title(f'grad($m$) against $\Delta$ ($m = n(\Delta-b)$): $(n, \
            b)$ = {fittings}', useTex = True)

    ax.set_ylabel(r'$m$', useTex = True)
    ax.set_xlabel(r'$\Delta$', useTex = True)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.grid(True, which='major', axis='both')
    fig.savefig(os.path.join(DIRECTORY_WRITE, 'optimal_freq.png'), format= 'png')
    plt.show(block=True)
    

if __name__ == '__main__':
    main()