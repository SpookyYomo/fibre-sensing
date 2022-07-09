# Gets the amplitude and angle deviation of a signal.

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
from scipy.signal import periodogram
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

# periodogram specific function
# specify sharing axes (https://stackoverflow.com/a/51684195)
def set_share_axes(axs, target=None, sharex=False, sharey=False):
    if target is None:
        target = axs.flat[0]
    # Manage share using grouper objects
    for ax in axs.flat:
        if sharex:
            target._shared_axes['x'].join(target, ax)
        if sharey:
            target._shared_axes['y'].join(target, ax)
    # Turn off x tick labels and offset text for all but the bottom row
    if sharex and axs.ndim > 1:
        for ax in axs[:-1,:].flat:
            ax.xaxis.set_tick_params(which='both', labelbottom=False, labeltop=False)
            ax.xaxis.offsetText.set_visible(False)
    # Turn off y tick labels and offset text for all but the left most column
    if sharey and axs.ndim > 1:
        for ax in axs[:,1:].flat:
            ax.yaxis.set_tick_params(which='both', labelleft=False, labelright=False)
            ax.yaxis.offsetText.set_visible(False)

# Place script specific functions here
def my_func():
    return

# For each file in files selected
def per_file(file, wd, display_plot, **kwargs):
    print() # New line
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
    
    # Get phases, amplitudes over time
    meta, trace = fr.parse_and_read_oscilliscope_txt(file)
    signal = signal_from_trace(np.asarray(trace))
    phases = signal_to_phase(signal, N, ph_ad, phase_advancement_correction= False)
    phases = phase_reconstruction_2(phases, ph_ad)
    amplitudes = get_R_signal(signal, N, ph_ad)
    t_axis = np.arange(start= 0, 
        stop= (int(meta["Record Length"][0])-N+1) * meta['Sample Interval'][0], \
        step= meta['Sample Interval'][0])
    time_axis = timeaxis_from_trace(trace)

    # Draw figure
    fig, axs = plt.subplots(nrows= 2, ncols= 3, sharex = False)
    set_share_axes(axs[0,:], sharex=True)
    set_share_axes(axs[1, 1:], sharex=True)
    
    # signal plot
    ax = axs[0, 0]
    try:
        ax.plot(time_axis, signal, color = 'mediumblue')
    except ValueError:
        dt = meta['Sample Interval'][0]
        axisT = np.fromiter([i* dt for i in range(len(signal))], 
            dtype= np.float64, count= len(signal))
        ax.plot(axisT, signal, color = 'mediumblue')
    ax.set_ylabel(r'$y(t)$/V', usetex= True)
    ax.set_xlabel(r'$t$/s', usetex= True)
    ax.set_title(f"Signal", usetex= False)

    # phase plot
    ax = axs[0, 1]
    ax.plot(t_axis, phases, color = 'mediumblue')
    # ax.plot(t_axis, linear(t_axis, *popt), color = 'darkgreen', markersize = 3)
    ax.set_ylabel(r'$\phi$/rad', usetex= True)
    ax.set_xlabel(r'$t$/s', usetex= True)
    ax.set_title(f"Phase", usetex= False)

    # amplitudes plot
    ax = axs[0, 2]
    ax.plot(t_axis, amplitudes, color = 'mediumblue')
    ax.set_ylabel(r'$A(t)$/V', usetex= True)
    ax.set_xlabel(r'$t$/s', usetex= True)
    ax.set_title(f"Amplitude", usetex= False)
    
    # spectral of signal plot
    ax = axs[1, 0]
    f, power_den = periodogram(signal, SAMPLING_F)
    print('\nSpectral Signal')
    print(f"Largest/peak frequency found at {f[np.argmax(power_den)]/1e3} kHz")
    ax.semilogy(f[0:]/1e3, power_den[0:])
    # ax.set_xscale('log')
    ax.set_ylim(min(power_den[1:])/5, max(power_den[1:])*5)
    ax.set_xlabel('frequency [kHz]')
    ax.set_ylabel('PSD [V**2/Hz]')
    ax.set_title('Signal FFT')

    # spectral of phase plot
    ax = axs[1, 1]
    f, Qxx_den = periodogram(phases, SAMPLING_F)
    print('\nSpectral Phase')
    print(f"Largest/peak frequency found at {f[np.argmax(Qxx_den)]/1e3} kHz")
    ax.semilogy(f[0:]/1e3, Qxx_den[0:])
    # ax.set_xscale('log')
    ax.set_ylim(min(Qxx_den[1:])/5, max(Qxx_den[1:])*5)
    ax.set_xlabel('frequency [kHz]')
    ax.set_ylabel('PSD [rad**2/Hz]')
    ax.set_title(f'phase FFT\nPeak at {f[np.argmax(Qxx_den)]/1e3 :.6f} kHz')
    
    # spectral of amplitude plot
    ax = axs[1, 2]
    f, Pxx_den = periodogram(amplitudes, SAMPLING_F)
    print('Spectral Amplitude')
    print(f"Largest/peak frequency found at {f[np.argmax(Pxx_den)]/1e3} kHz")
    ax.semilogy(f[1:]/1e3, Pxx_den[1:])
    # ax.set_xscale('log')
    ax.set_ylim(min(Pxx_den[1:])/5, max(Pxx_den[1:])*5)
    ax.set_xlabel('frequency [kHz]')
    ax.set_ylabel('PSD [V**2/Hz]')
    ax.set_title('amplitudes FFT')

    # zoom in to observe 200Hz amplitude oscillations
    axs[0, 2].set_xlim(0, 2/200)

    fig.suptitle(f"Sampling freq $f= ${SAMPLING_F/1e6} MHz, AOM freq = \
        {SIGNAL_F/1e6:.6f}MHz, N = {N}", usetex= True)
    # powerpoint is 13.333 inches wide by 7.5 inches high
    # a4 paper is 8.25 by 11.75 in
    # geometry package boundary -0.875in*2 horizontally
    # figures set to 0.8 * textwidth
    # fig.set_size_inches(0.8*(8.25-0.875*2), 2) 
    fig.set_size_inches(11.75-1.5, 8.25-2 - 1)
    fig.tight_layout()

    if display_plot:
        plt.show(block= True)
        plt.pause(0.8)

    # Saving figure block
    try:
        print('Saving plot...')
        fig.savefig(os.path.join(wd, NAME[:-4] + 'result.png'), 
            format= 'png', dpi= 300)
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
        per_file(file, DIR_WRITE, len(files) == 1)
    print("Script has ended.")

if __name__ == '__main__':
    main()