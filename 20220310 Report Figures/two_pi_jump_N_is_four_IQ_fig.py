# Generate 2 pi jump when N = 4 from 20211214 Linearity with Amp and Freq Incomplete\Batch1\C1-004Hz-020.4mVpp-80.125000MHz-OneAmp00001.txt

# . Honours Module Folder
# ├ FYPLibrary
# | ├ file_reading.py
# | └ IQ_demod.py
# └ Sub Project folder
#   ├ traces folder
#   | └ C1-004Hz-020.4mVpp-80.125000MHz-OneAmp00001.txt
#   └ this script(.py)

# Initialiastion: Directory appending for my system. Vary the directories as necessary.
import sys, os.path

# Add the FYPLibrary folder in the same level as this folder.
if os.path.dirname(os.path.dirname(__file__)) not in sys.path:
    sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'FYPLibrary'))

# from pprint import pprint
# pprint(sys.path)

# Import Modules
from os import listdir
from file_reading import *
from IQ_demod import *
import numpy as np
# from numpy import pi
# import EPstandard
# from scipy.signal import periodogram
import matplotlib
import matplotlib.pyplot as plt
# from matplotlib.ticker import AutoMinorLocator
import string
from scipy.interpolate import interp1d

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
    
def polar(x, y) -> tuple:
    """returns rho, theta (degrees)"""
    return np.hypot(x, y), np.arctan2(y, x)

def main():
    file = os.path.join(os.path.dirname(__file__), 
        "traces", "C1-004Hz-020.4mVpp-80.125000MHz-OneAmp00001.txt")
    print(f"File selected: {file}")
    wd = write_directory(__file__)
    print(f"Results will be written to: {wd}")

    # Edit Initialisation
    SIGNAL_F = 80.125e6*2 #Hz 
    SAMPLING_F = 1.0e6 #Hz
    ph_ad = phase_advance(SIGNAL_F, SAMPLING_F) # phase advance = 2*pi/N
    N, _ = freq_ratio(signal=SIGNAL_F, sample=SAMPLING_F)
    print(f"[Int Debug] {N = }, {ph_ad*2/np.pi = } should be 1")

    meta, trace = parse_and_read_oscilliscope_txt(file)
    signal = signal_from_trace(trace)
    phases = signal_to_phase(signal, N, ph_ad, phase_advancement_correction= False)
    phases = phase_reconstruction_2(phases, ph_ad)
    t_axis = np.arange(start= 0, 
        stop= (int(meta["Record Length"][0])-N+1) * meta['Sample Interval'][0], step= meta['Sample Interval'][0])
    # radii = get_R_signal(signal, N, ph_ad)

    zoom_in_start = 113695
    zoom_in_end = zoom_in_start + 7

    matplotlib.rcParams.update({'font.size': 11})

    fig = plt.figure()
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 2])

    # Voltage over time
    ax0 = fig.add_subplot(gs[0, 0])
    ax0.plot(t_axis[zoom_in_start:zoom_in_end+6]*1e6, 
        signal[zoom_in_start:zoom_in_end+6]*1e3, color = 'mediumblue')
    ax0.set_ylabel(r'$V$/mV', usetex= True)
    ax0.set_ylim([min(phases) - 0.5, max(phases) + 0.5])
    # ax0.set_ylim(min(signal[zoom_in_start:zoom_in_end+6]), max(signal[zoom_in_start:zoom_in_end+6])+)
    ax0.set_ylim(np.asarray([-0.05, 0.05])*1e3)
    ax0.axes.xaxis.set_visible(False)

    # phase slip
    ax1 = fig.add_subplot(gs[1, 0])
    ax1.plot(t_axis[zoom_in_start:zoom_in_end+6]*1e6-t_axis[zoom_in_start]*1e6, 
        phases[zoom_in_start:zoom_in_end+6], color = 'mediumblue')
    ax1.set_ylabel(r'$\phi_d$/rad', usetex= True)
    ax1.set_xlabel(r'$t$/$\mu$s', usetex= True)
    # ax1.axes.xaxis.set_visible(False)

    # circle plot
    Is, Qs = get_Is_Qs(signal, N, ph_ad)
    Is = Is[zoom_in_start+N:zoom_in_end+N]
    Qs = Qs[zoom_in_start+N:zoom_in_end+N]
    Rs, thetas = polar(Is, Qs)
    ax2 = fig.add_subplot(gs[:, 1], projection= 'polar')
    ax2.scatter(thetas, Rs, c= 100*Rs/np.max(Rs), cmap= 'coolwarm')
    for i, t in enumerate(zip(thetas, Rs)):
        ax2.annotate(i+1, t)
    ax2.set_rticks([0.010, 0.02, 0.03], minor= False)
    # ax.set_theta_offset(-np.pi)
    # ax.set_thetalim(-np.pi / 1, np.pi / 1.5)
    xT=plt.xticks()[0]
    xL=['0',r'$\frac{\pi}{4}$',r'$\frac{\pi}{2}$',r'$\frac{3\pi}{4}$',\
    r'$\pi$',r'$\frac{5\pi}{4}$',r'$\frac{3\pi}{2}$',r'$\frac{7\pi}{4}$']
    plt.xticks(xT, xL)

    my_zip = list(zip(thetas, Rs))
    for i in range(len(my_zip)-1):
        (th1, r1) = my_zip[i]
        (th2, r2) = my_zip[i+1]
        if th2 < th1: th2 += 2*np.pi
        th_space = np.linspace( th1, th2, 100)
        r_space = np.linspace( r1, r2, 100)
        ax2.plot(th_space, r_space, color='black', alpha= 0.3 )

    axs = [ax0, ax1, ax2]
    for n, ax in enumerate(axs):  
        ax.text(-0.1, 1.1, "(" + string.ascii_lowercase[n] + ")", transform=ax.transAxes, 
                size=11, weight='bold')

    # powerpoint is 13.333 inches wide by 7.5 inches high
    # geometry package boundary -0.875in*2 horizontally
    # figures set to 0.85 * textwidth
    fig.set_size_inches(0.85*(8.25-0.875*2), 3.60) 
    fig.tight_layout()

    # plt.show(block= True)
    # # Saving figure block
    try:
        print('Saving plot...')
        NAME = os.path.basename(file)
        fig.savefig(os.path.join(wd, 'IQ_plane_' + NAME[:-4] + '_2pijumpIQ.png'), 
            format= 'png', dpi= 300, pad_inches = 0.01)
        print('Saved')
    except FileNotFoundError:
        print(f"[Error] Folder should have been created in initialisation block")
        exit()
    plt.close()

if __name__ == '__main__':
    main()