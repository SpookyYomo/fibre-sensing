# Initialiastion: Directory appending for my system. Vary the directories as necessary.
import sys
import os.path
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
import scipy
from scipy import optimize

# Yield all files in batch 4
files = listdir(path=os.path.dirname(__file__)+'\Batch 4')

# This file specific function
def get_f_and_A_from_filename(txt):
    Hz_pos = txt.find('Hz')
    mVpp_pos = txt.find('mVp')
    f = int(txt[Hz_pos-3:Hz_pos])
    A = int(txt[mVpp_pos-3:mVpp_pos])
    return f, A

def fit_fn(x, A, f, phi, m, c):
    return A * np.sin(2*pi*f*x + phi) + m*x+c

SIGNAL_F = 80.125e6*2 #Hz 
SAMPLING_F = 1.0e6 #Hz
N, M = freq_ratio(signal=SIGNAL_F, sample=SAMPLING_F)

results = []
for file in files:
    NAME = file
    TXT_FILE_PATH = file_str_to_path(NAME, dir=os.path.dirname(__file__)+'\Batch 4')
    
    f, mVpp = get_f_and_A_from_filename(NAME)
    if f == 0:
        continue
    print(f"{NAME = }, {f = }, {mVpp = }")

    readTxt = fr.read_oscilliscope_txt(TXT_FILE_PATH)
    meta, trace = fr.parse_oscilliscope_txt(readTxt)
    signal = signal_from_trace(np.asarray(trace))
    phases = signal_to_phase(signal, N, 2*pi/N)
    phases = phase_reconstruction(phases, 5)
    t_axis = np.arange(start= 0, 
        stop= (int(meta["Record Length"][0])-N) * meta['Sample Interval'][0], step= meta['Sample Interval'][0])
    
    popt, pcov = optimize.curve_fit(fit_fn, t_axis, phases, p0 = [2.5/50 * mVpp, f, 4.5, -80, 0])
    fittings = EPstandard.easy_read_popt_pcov(popt, pcov)
    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.plot(t_axis, phases, color = 'mediumblue')
    ax.plot(t_axis, fit_fn(t_axis, *popt), color = 'darkgreen', markersize = 3)
    ax.set_ylabel(r'$\phi$/rad', useTex = True)
    ax.set_xlabel(r'$t$/s', useTex = True)
    plt.title(f"Mirror Freq = {f} Hz, Piezo Ampl = {mVpp} mVpp, Signal Freq = {SIGNAL_F:6.0F}Hz, \n(A, $f$, $\phi_0$, m, c) = {fittings}", useTex = True )
    plt.show(block = False)
    plt.pause(0.4)
    DIRECTORY_WRITE = os.path.dirname(__file__)
    fig.savefig(os.path.join(DIRECTORY_WRITE, NAME[:-4] + '.png'))
    plt.close()
    
    results.append((mVpp, popt[0], np.sqrt(pcov[0][0])))

xs = [r[0] for r in results]
ys = [abs(r[1]) for r in results]
yerrs = [r[2] for r in results]

fig, ax = plt.subplots(nrows=1, ncols=1)
ax.errorbar(xs, ys, 
    yerr=yerrs, fmt = 'bo', markersize = 3)
ax.set_ylabel(r'$A$/rad', useTex = True)
ax.set_xlabel(r'mVpp', useTex = True)
plt.title('A against mVpp')
plt.show(block=True)