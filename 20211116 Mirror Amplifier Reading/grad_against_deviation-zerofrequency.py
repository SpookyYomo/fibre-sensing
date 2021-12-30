# Initialiastion: Directory appending for my system. Vary the directories as necessary.
import sys, os.path
# Add the FYPLibrary folder in the same level as this folder.
if os.path.dirname(os.path.dirname(__file__)) not in sys.path:
    # print(f"{os.path.dirname(os.path.dirname(__file__)) = }")
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'FYPLibrary'))

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

# Yield all files in batch folder
folderpath = os.path.dirname(__file__)+'\Batch 2b'
files = listdir(path=folderpath)

# This file's specific function
def linear(x, m, r):
    return m*(x-r)

def main():
    SIGNAL_F = 80.125e6*2 #Hz 
    SAMPLING_F = 1.0e6 #Hz
    N, M = freq_ratio(signal=SIGNAL_F, sample=SAMPLING_F)
    INDIVIDUAL_PLOTS = True
    DISPLAY_INDIVIDUAL_PLOTS = False
    FINAL_REGRESSION = True
    eighty_million_displaced = True # For recording delta_corr as absolute or relative to 80.125MHz

    results = []
    for file in files:
        NAME = file
        TXT_FILE_PATH = file_str_to_path(NAME, dir=folderpath)
        if NAME[0] != 'C':
            print(f"Unrecognized file type by name of {NAME}")
            continue
        
        # delta corr should have the unit of Hz here.
        # f, mVpp, delta_corr = get_f_A_and_deviation_from_filename(NAME)
        _, f, mVpp, delta_corr = numerical_variables_from_name(NAME)
        if f != 0:
            print('file with oscillating mirror was encountered!')
            continue
        if eighty_million_displaced:
            delta_corr -= 80.125
        delta_corr *= 1e6 # converts MHz to Hz
        print(f"[Int Debug]{NAME = }, {f = }, {mVpp = }, {delta_corr = }")
        # print(f"[Int Debug]{float('80.1250'+str(delta_corr)) = }")

        # readTxt = fr.read_oscilliscope_txt(TXT_FILE_PATH)
        # meta, trace = fr.parse_oscilliscope_txt(readTxt)
        meta, trace = fr.parse_and_read_oscilliscope_txt(TXT_FILE_PATH)
        signal = signal_from_trace(np.asarray(trace))
        phases = signal_to_phase(signal, N, 2*pi/N, True)
        phases = phase_reconstruction(phases, 4.8)
        t_axis = np.arange(start= 0, 
            stop= (int(meta["Record Length"][0])-N) * meta['Sample Interval'][0], step= meta['Sample Interval'][0])
        
        # set initial to 0 value
        phases -= phases[0]
        popt, pcov = optimize.curve_fit(linear, t_axis, phases, p0=[10*delta_corr, 0] )
        fittings = EPstandard.easy_read_popt_pcov(popt, pcov)
        
        if INDIVIDUAL_PLOTS:
            fig, ax = plt.subplots(nrows=1, ncols=1)
            ax.plot(t_axis, phases, color = 'mediumblue')
            ax.plot(t_axis, linear(t_axis, *popt), color = 'darkgreen', markersize = 3)
            ax.set_ylabel(r'$\phi$/rad', useTex = True)
            ax.set_xlabel(r'$t$/s', useTex = True)
            plt.title(f"$\Delta$ = {delta_corr:.5f}Hz, ($m$, $t_0$) = {fittings}", useTex = True )
            if DISPLAY_INDIVIDUAL_PLOTS:
                plt.show(block = False)
                plt.pause(0.8)
            DIRECTORY_WRITE = os.path.dirname(__file__)
            fig.savefig(os.path.join(DIRECTORY_WRITE, NAME[:-4] + '.png'))
            plt.close()
        
        results.append((delta_corr, popt[0], np.sqrt(pcov[0][0])))

    xs = [r[0] for r in results]
    ys = [r[1] for r in results]
    yerrs = [r[2] for r in results]

    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.errorbar(xs, ys, 
        yerr=yerrs, fmt = 'bo', markersize = 3)

    if FINAL_REGRESSION:
        popt, pcov = optimize.curve_fit(linear, xs, ys, p0=(100, 55), sigma= yerrs)
        fittings = EPstandard.easy_read_popt_pcov(popt,pcov)
        ax.plot(delta_axis:= np.arange(start=np.min(xs), stop=np.max(xs), step = 1), linear(delta_axis, *popt),
            markersize = 1)
        plt.title(f'grad($m$) against $\Delta$ ($m = n(\Delta-b)$): $(n, b)$ = {fittings}', useTex = True)

    ax.set_ylabel(r'$m$', useTex = True)
    ax.set_xlabel(r'$\Delta$', useTex = True)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.grid(True, which='major', axis='both')
    
    plt.show(block=True)

if __name__ == '__main__':
    main()