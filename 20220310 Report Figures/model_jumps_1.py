# Generate toy model

# . Honours Module Folder
# ├ FYPLibrary
# | ├ file_reading.py
# | └ IQ_demod.py
# ├ Report Figures
# | └ this script(.py)
# └ 20220316 130 meter cable
#   ├ Batch 2
#   . └ C1-1.0MHz-00000.trc


# Initialiastion: Directory appending for my system. Vary the directories as necessary.
import sys, os.path

# Add the FYPLibrary folder in the same level as this folder.
if os.path.dirname(os.path.dirname(__file__)) not in sys.path:
    sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'FYPLibrary'))

# Import Modules
from file_reading import *
from IQ_demod import *
import numpy as np
from numpy import random as rd
from numpy import pi
# from numpy import pi
# import EPstandard
# from scipy.signal import periodogram
import matplotlib
import matplotlib.pyplot as plt
# from matplotlib.ticker import AutoMinorLocator
from scipy.signal import periodogram
import string

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

def main():
    ####### Create toy model 2pi jumps
    # file = os.path.join(os.path.dirname(os.path.dirname(__file__)), 
    #     "20220316 130 meter cable", "Batch 4" ,"Batch 4F 8", "C1-1.0MHz-00001.trc")
    # print(f"File selected: {file}")
    wd = os.path.join(write_directory(__file__), "Model_Jumps")
    print(f"Results will be written to: {wd}")

    # init
    rd.seed(1337) # seed random
    SIGNAL_F = 80.125e6*2 #Hz 
    SAMPLING_F = 1.0e6 #Hz
    N, _ = freq_ratio(signal=SIGNAL_F, sample=SAMPLING_F)
    ph_ad = phase_advance(SIGNAL_F, SAMPLING_F) # phase advance = 2*pi/N
    print(f"[Int Debug] {N = }, {ph_ad = }")
    amp_V = 0.092 # amplitude of V
    phi_0 = 0.163 # radian associated with t = 0
    mirror_f = 0.5 # freq of oscillation of mirror
    mirror_A = 0.8 # amplitude of oscillating mirror
    rel_noise_amp = 0.02 # amplitude of noise of detector 
    rel_noise_std = 0.3 # std deviation
    phase_noise_amp = 1 # noise of mirror amplitude relative to miror_A
    phase_std = 0.48 # std_deviation
    t_max = 4 #s

    t_axis = np.arange(0, t_max, 1/SAMPLING_F)
    len_t_axis = len(t_axis)

    # model 2
    signal = amp_V * np.sin(2*pi*SIGNAL_F * t_axis + phi_0 \
        + mirror_A * (np.sin(2*pi*mirror_f * t_axis) + \
            phase_noise_amp * rd.default_rng(13).normal(0, phase_std, size= len_t_axis))) + \
        rel_noise_amp * rd.default_rng(2).normal(0, rel_noise_std, size= len_t_axis)

    phases = signal_to_phase(signal, N, ph_ad, phase_advancement_correction= False)
    phases = phase_reconstruction_2(phases, ph_ad)

    matplotlib.rcParams.update({'font.size': 11})

    fig, ax = plt.subplots(nrows=1, ncols=1, sharex= True)
    ax.plot(t_axis[N-1:], phases, color = 'mediumblue', linewidth= 0.1)
    ax.set_ylabel(r'$\phi_d$/rad', usetex= True)
    ax.set_xlabel(r'$t$/s', usetex= True)
    ax.set_ylim(-1.5*2*pi, 1.5*2*pi)
    ax.set_xlim([0.0, 4])

    fig.set_size_inches(0.85*(8.25-0.875*2), 3) 
    fig.tight_layout()

    # Saving figure block
    try:
        print('Saving plot...')
        fig.savefig(os.path.join(wd, "model_jump_N_4.png"), 
            format= 'png', dpi= 300, pad_inches = 0.00)
        print('Saved')
    except FileNotFoundError:
        print(f"[Error] Folder should have been created in initialisation block")
        exit()
    plt.close()

    ###### Power spectrum density (toy model vs real data)
    
    # toy model
    f1, power_den1 = periodogram(signal, SAMPLING_F)

    # plot toy model
    fig, axs = plt.subplots(nrows= 1, ncols= 2, sharey = True)
    ax = axs[0]
    ax.semilogy(f1[0:]/1e3, power_den1[0:], linewidth= 0.1)
    # ax.set_xscale('log')
    # ax.set_ylim(min(power_den[1:])/5, max(power_den[1:])*5)
    ax.set_xlabel(r'$f$ /kHz', usetex= True)
    ax.set_ylabel(r'PSD /V${ }^2$ Hz${ }^{-1}$', usetex= True)

    # obtain real data and process
    file = os.path.join(os.path.dirname(__file__), 
        "traces", "C1-004Hz-020.4mVpp-80.125000MHz-OneAmp00001.txt")
    print(f"File selected: {file}")
    SIGNAL_F = 80.125e6*2 #Hz 
    SAMPLING_F = 1.0e6 #Hz
    ph_ad = phase_advance(SIGNAL_F, SAMPLING_F) # phase advance = 2*pi/N
    N, _ = freq_ratio(signal=SIGNAL_F, sample=SAMPLING_F)
    _, trace = fr.parse_and_read_oscilliscope_txt(file)
    signal = signal_from_trace(np.asarray(trace))

    # plot real data
    f2, power_den2 = periodogram(signal, SAMPLING_F)
    ax = axs[1]
    ax.semilogy(f2[0:]/1e3, power_den2[0:], linewidth= 0.1)
    # ax.set_xscale('log')
    ax.set_ylim(min(min(power_den1[1:]), min(power_den2[1:]))/5, 
        max(max(power_den1[1:]), max(power_den2[1:]))*5)
    ax.set_xlabel(r'$f$ /kHz', usetex= True)
    ax.set_ylabel(r'PSD /V${ }^2$ Hz${ }^{-1}$', usetex= True)
    ax.axes.yaxis.set_visible(False)

    # format subplots
    for n, ax in enumerate(axs.flat):  
        ax.text(-0.1, 1.1, "(" + string.ascii_lowercase[n] + ")", transform=ax.transAxes, 
                size=11, weight='bold')
    fig.set_size_inches(0.85*(8.25-0.875*2), 3) 
    fig.tight_layout()
    try:
        print('Saving plot...')
        fig.savefig(os.path.join(wd, "model_vs_real_PSD_N_4.png"), 
            format= 'png', dpi= 300, pad_inches = 0)
        print('Saved')
    except FileNotFoundError:
        print(f"[Error] Folder should have been created in initialisation block")
        exit()
    plt.close()


    ###### Create for N = 8
    # init
    rd.seed(1337) # seed random
    SIGNAL_F = 80.0625e6*2 #Hz 
    SAMPLING_F = 1.0e6 #Hz
    N, _ = freq_ratio(signal=SIGNAL_F, sample=SAMPLING_F)
    ph_ad = phase_advance(SIGNAL_F, SAMPLING_F) # phase advance = 2*pi/N
    print(f"[Int Debug] {N = }, {ph_ad = }")
    amp_V = 0.092 # amplitude of V
    phi_0 = 0.163 # radian associated with t = 0
    mirror_f = 0.5 # freq of oscillation of mirror
    mirror_A = 0.8 # amplitude of oscillating mirror
    rel_noise_amp = 0.02 # amplitude of noise of detector 
    rel_noise_std = 0.3 # std deviation
    phase_noise_amp = 1 # noise of mirror amplitude relative to miror_A
    phase_std = 0.48 # std_deviation
    t_max = 4 #s

    t_axis = np.arange(0, t_max, 1/SAMPLING_F)
    len_t_axis = len(t_axis)

    # model 2
    signal = amp_V * np.sin(2*pi*SIGNAL_F * t_axis + phi_0 \
        + mirror_A * (np.sin(2*pi*mirror_f * t_axis) + \
            phase_noise_amp * rd.default_rng(13).normal(0, phase_std, size= len_t_axis))) + \
        rel_noise_amp * rd.default_rng(2).normal(0, rel_noise_std, size= len_t_axis)

    phases = signal_to_phase(signal, N, ph_ad, phase_advancement_correction= False)
    phases = phase_reconstruction_2(phases, ph_ad)

    fig, ax = plt.subplots(nrows=1, ncols=1, sharex= True)
    ax.plot(t_axis[N-1:], phases, color = 'mediumblue', linewidth= 0.1)
    ax.set_ylabel(r'$\phi_d$/rad', usetex= True)
    ax.set_xlabel(r'$t$/s', usetex= True)
    ax.set_ylim(-3, 3)
    ax.set_xlim([0.0, 4])

    fig.set_size_inches(0.85*(8.25-0.875*2), 3) 
    fig.tight_layout()

    # Saving figure block
    try:
        print('Saving plot...')
        fig.savefig(os.path.join(wd, "model_jump_N_8.png"), 
            format= 'png', dpi= 300, pad_inches = 0.00)
        print('Saved')
    except FileNotFoundError:
        print(f"[Error] Folder should have been created in initialisation block")
        exit()
    plt.close()

    return

if __name__ == '__main__':
    main()