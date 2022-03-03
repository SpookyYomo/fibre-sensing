# When the oscilliscope is not using a reference clock, tweaking the AOM 
# frequency is a viable means of minimising the extent of drift. This 
# script aims to help pin down the appropriate deviation to set the AOM.

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
import EPstandard
from scipy import optimize
# from scipy.signal import periodogram
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

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
    my_input = input(f"Results will be written to: {wd}\nAccept? [Y/N]")
    if my_input == 'N':
        import tkinter as tk
        from tkinter import filedialog
        wd = filedialog.askdirectory()
        print(f"Results will be written to {wd}")
    if not os.path.exists(wd):
        os.makedirs(wd)
    return wd

# Script specific functions here
def linear(x, m, c):
    return m*x+c

def drifting_sin(x, A, f, phi, m, c):
    return A * np.sin(2*pi*f*x + phi) + m*x+c

# For each file in files selected
def per_file(file, wd, gen_plot, display_plot):
    # wd is the write directory for results of this given file
    # Edit kwargs as necessary
    
    # Read and obtain trace
    NAME = os.path.basename(file)
    if NAME[0] != 'C' and NAME[-4:] != '.txt':
        print(f"[Warning] Unrecognised file reached, Skipping: {NAME}\n")
        return

    # Initialisation
    _, f, mVpp, aom_freq = numerical_variables_from_name(NAME)
    SIGNAL_F = 80.125e6*2 #Hz 
    SAMPLING_F = 1.0e6 #Hz
    ph_ad = phase_advance(SIGNAL_F, SAMPLING_F) # phase advance = 2*pi/N
    N, _ = freq_ratio(signal= SIGNAL_F, sample= SAMPLING_F)
    delta_corr = (aom_freq-80.125)*1e6 # gets Hz difference from 80.125MHz
    print(f"[Int Debug] {NAME = }, {f = }, {mVpp = }, delta_corr = {delta_corr:.3f}Hz")
    
    # Get phases over time
    meta, trace = fr.parse_and_read_oscilliscope_txt(file)
    signal = signal_from_trace(np.asarray(trace))
    phases = signal_to_phase(signal, N, ph_ad, phase_advancement_correction= False)
    phases = phase_reconstruction_2(phases, ph_ad)
    t_axis = np.arange(start= 0, 
        stop= (int(meta["Record Length"][0])-N+1) * meta['Sample Interval'][0], step= meta['Sample Interval'][0])

    # perform fitting
    drifting_sin_bounds = ((     0,      0,  -2*pi, -np.inf, -np.inf), \
                           (np.inf, np.inf,   2*pi,  np.inf,  np.inf))
    popt, pcov = optimize.curve_fit(drifting_sin, t_axis, phases, 
        p0= [0.05*mVpp, f, 0, delta_corr*12.4, 0.5*(max(phases[:1000])+min(phases[:1000]))], 
        bounds= drifting_sin_bounds)
    fittings = EPstandard.easy_read_popt_pcov(popt, pcov)

    # Draw figure
    if gen_plot:
        fig, ax = plt.subplots(nrows=1, ncols=1)
        ax.plot(t_axis, phases, color = 'mediumblue')
        ax.set_ylabel(r'$\phi$/rad', usetex= True)
        ax.set_xlabel(r'$t$/s', usetex= True)
        plt.title(f"$f_{{Mirror}}$ = {f}Hz, Piezo Ampl = {mVpp}mVpp, \
            Signal Freq = {SIGNAL_F:6.0F}Hz, \n(A, $f$, $\phi_0$, m, \
            c) = {fittings}", usetex= True)
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

    return (mVpp, popt[0], np.sqrt(pcov[0][0]))

def final_movement(accum, dw):
    # perform a of linearity against deviation
    xs = [r[0] for r in accum]
    ys = [abs(r[1]) for r in accum]
    yerrs = [r[2] for r in accum]
    print("Commencing final diagram...")
    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.errorbar(xs, ys, 
        yerr=yerrs, fmt = 'bo', markersize = 3)
    
    lin_bounds = ((0, -np.inf), (np.inf, np.inf))
    popt, pcov = optimize.curve_fit(linear, xs, ys, p0=(0.05, 0), 
        sigma= yerrs, bounds= lin_bounds)
    fittings = EPstandard.easy_read_popt_pcov(popt, pcov)
    ax.plot(mVpp_axis:= np.arange(start=np.min(xs), stop=np.max(xs), \
            step = 0.005*(np.max(xs)-np.min(xs))), 
        linear(mVpp_axis, *popt),
        markersize = 1)
    plt.title(f'$A$ against mVpp\n ($A$ = $m$*mVpp + $c$: $(m, c)$ = \
        {fittings}', usetex= True)

    ax.set_ylabel(r'$A$/rad', usetex= True)
    ax.set_xlabel(r'mVpp', usetex= True)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.grid(True, which='major', axis='both')
    fig.savefig(os.path.join(dw, 'linearity.png'), format= 'png')
    plt.show(block=True)
    return

def main():
    print("Select Files to perform this script on.")
    files = get_files()
    print(f"Files selected: {files}")
    DIR_WRITE = write_directory(files)
    print(f"Results will be written to: {DIR_WRITE}")
    results = []
    for file in files:
        results.append(per_file(file, DIR_WRITE, True, False))
    final_movement(results, DIR_WRITE)
    print("Script has ended.")

if __name__ == '__main__':
    main()