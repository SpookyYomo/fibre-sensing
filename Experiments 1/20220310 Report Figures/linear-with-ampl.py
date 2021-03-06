# This script will performt the regression over all traces opened, and 
# perform a drifting sinusoidal fit. All the collected results are then 
# bunched together to see if they obey a linear relationship with 
# mirror_amplitude. The figure that we generate takes trace data from

# . Honours Module Folder
# ├ FYPLibrary
# | ├ file_reading.py
# | └ IQ_demod.py
# └ Experiments 1 (Optional Depth)
#   └ Sub Project folder
#     ├ Batch Folders
#     | └ C1-...-000001.txt
#     └ this script(.py)

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
import EPstandard
from scipy import optimize
# from scipy.signal import periodogram
import matplotlib.pyplot as plt
# from matplotlib.ticker import AutoMinorLocator
import matplotlib
import string

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
def linear(x, m, c):
    return m*x+c

def drifting_sin(x, A, f, phi, m, c):
    return A * np.sin(2*pi*f*x + phi) + m*x+c

# For each file in files selected
def per_file(enum_index, file, wd, gen_plot, display_plot, **kwargs):
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
    print(f"[Int Debug] {num_vars = }, {N = }")
    _, mirror_f, mVpp  = num_vars
    if mirror_f == 0:
        print("Skipping 0 frequency")
        return
    mVpp *= 50 # dds step up before entering into pieze

    # Get phases over time
    meta, trace = fr.parse_and_read_oscilliscope_txt(file)
    signal = signal_from_trace(np.asarray(trace))
    phases = signal_to_phase(signal, N, ph_ad, phase_advancement_correction= False)
    phases = phase_reconstruction_2(phases, ph_ad)
    t_axis = np.arange(start= 0, 
        stop= (int(meta["Record Length"][0])-N+1) * meta['Sample Interval'][0], step= meta['Sample Interval'][0])

    # Regression to obtain mirror freq
    drifting_sin_bounds = ((     0,      0,  -2*pi, -np.inf, -np.inf), \
                           (np.inf, np.inf,   2*pi,  np.inf,  np.inf))
    popt, pcov = optimize.curve_fit(drifting_sin, t_axis, phases, 
        p0= [0.05*mVpp/50, mirror_f, 0, 0, 0.5*(max(phases[:1000])+min(phases[:1000]))], 
        bounds= drifting_sin_bounds)
    # fittings = EPstandard.easy_read_popt_pcov(popt, pcov)

    # Subtract drift for figure
    phases -= popt[-2]*t_axis

    # # Regression to obtain mirror freq
    # drifting_sin_bounds = ((     0,      0,  -2*pi, -np.inf, -np.inf), \
    #                        (np.inf, np.inf,   2*pi,  np.inf,  np.inf))
    # popt, pcov = optimize.curve_fit(drifting_sin, t_axis, phases, 
    #     p0= [0.05*mVpp/50, mirror_f, 0, 0, 0.5*(max(phases[:1000])+min(phases[:1000]))], 
    #     bounds= drifting_sin_bounds)
    fittings = EPstandard.easy_read_popt_pcov(popt, pcov)
    
    # we pass this value into main() for appending
    value = (mVpp, popt[0], np.sqrt(pcov[0][0])) 

    # Draw figure
    # if False:
    if enum_index == 3 or enum_index == 8:
    # if enum_index:
        fig, ax = plt.subplots(nrows=1, ncols=1)
        ax.plot(t_axis, phases, color = 'mediumblue', linewidth= 0.1)
        ax.set_ylabel(r'$\phi_d$/rad', usetex= True)
        ax.set_xlabel(r'$t$/s', usetex= True)
        ax.set_xlim([0, 0.2])
        ax.text(-0.1, 1.1, "(" + string.ascii_lowercase[0] + ")", transform=ax.transAxes, 
                size=11, weight='bold')
        # plt.title(f"$f_{{Mirror}}$ = {mirror_f} Hz, Piezo Ampl = {mVpp/1e3} Vpp, \
        #     \n $f_s$ = {SAMPLING_F/1e6:.2f} MS/s, $2\Omega$ = {SIGNAL_F/1e6:.3f}MHz, $N$ = {N} \
        #     \n(A, $f$) = {fittings[:2]}",
        #     usetex= True )
        # powerpoint is 13.333 inches wide by 7.5 inches high
        # geometry package boundary -1.0in*2 horizontally
        # figures set to 0.85 * textwidth
        # fig.set_size_inches(0.85*(8.25-1.0*2), 2.7) 
        fig.set_size_inches(4, 2.25) 
        fig.tight_layout()
        if display_plot:
            plt.show(block= True)
            plt.pause(0.8)
        # Saving figure block
        try:
            print('Saving plot...')
            fig.savefig(os.path.join(wd, 'oscillation_eg_' + NAME[:-4] + 'result_notitle.png'), 
                format= 'png', dpi = 300, pad_inches = 0)
            print('Saved')
        except FileNotFoundError:
            print(f"[Error] Folder should have been created in initialisation block")
            exit()
        plt.close()

    return value

def final_movement(result, wd):
    xs = [r[0] for r in result]
    ys = [abs(r[1]) for r in result]
    yerrs = [r[2] for r in result]

    xs=np.asarray(xs)/1e3 #mV to V 

    # Regress obtained values of fitted frequency against input frequency
    lin_bounds = ((0, -np.inf), (np.inf, np.inf))
    popt, pcov = optimize.curve_fit(linear, xs, ys, p0=(1, 0), 
        sigma= yerrs, bounds= lin_bounds)
    fittings = EPstandard.easy_read_popt_pcov(popt,pcov)

    # Plot out the regressed data
    fig, ax = plt.subplots(nrows=1, ncols=1)
    # plt.grid(True, which='major', axis='both')
    ax.errorbar(xs, ys, yerr=yerrs, fmt = 'bo', markersize = 3)
    ax.plot(xAx:= np.arange(start=np.min(xs), stop=np.max(xs), \
        step = 0.005*(np.max(xs)-np.min(xs))),
        linear(xAx, *popt),
        markersize = 1)
    # plt.title(f"Amplitude against mirror oscillation\n \
    #     $y=mx+c$: ($m, c$) = {fittings}")
    ax.set_ylabel(r'$A_{Regressed}$/rad', usetex= True)
    ax.set_xlabel(r'Vpp/V', usetex= True)
    ax.set_xlim([0, 20])
    ax.set_ylim([0, 20])
    ax.text(-0.1, 1.1, "(" + string.ascii_lowercase[1] + ")", transform=ax.transAxes, 
                size=11, weight='bold')

    # Save plot
    fig.set_size_inches(2.1, 2.25) 
    fig.tight_layout()
    try:
        print('Saving plot...')
        fig.savefig(os.path.join(wd, 'oscillation_linear_with_amplitude_notitle.png'), 
            format= 'png', dpi= 300, pad_inches = 0)
        print('Saved')
    except FileNotFoundError:
        print(f"[Error] Folder should have been created in initialisation block")
        exit()
    return

def main():
    print("Select Files to perform this script on.")
    hons_folder = os.path.dirname(os.path.dirname(__file__))
    files_folder = os.path.join(hons_folder, 
        "20211020 Linear With Amplt and Deviation", "Batch 4")
    files = [os.path.join(files_folder,f) for f in os.listdir(files_folder) if os.path.isfile(os.path.join(files_folder,f))]  
    print(f"Files selected: {files}")
    DIR_WRITE = write_directory(__file__)
    print(f"Results will be written to: {DIR_WRITE}")
    results = []
    matplotlib.rcParams.update({'font.size': 11})
    for enum_index, file in enumerate(files):
        temp_append = per_file(enum_index, file, DIR_WRITE, True, len(files)==1)
        if temp_append != None: results.append(temp_append)
    if len(set(map(lambda x: x[0], results))) <= 1: 
        print("Insufficient x-axis to determine linear relationship.\n \
            Script has ended")
        return
    print("\nWith all data processed, now performing final regression.")
    final_movement(results, DIR_WRITE)
    print("Script has ended.")

if __name__ == '__main__':
    main()