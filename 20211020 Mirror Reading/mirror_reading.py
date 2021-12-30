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
from file_reading import *
from IQ_demod import *
import numpy as np
import matplotlib.pyplot as plt

## Main code
def main():
    # choose the appropriate data file
    ## Batch 1
    # SIGNAL_F = 80.125e6*2 #Hz 
    # SAMPLING_F = 1.0e6 #Hz
    # NAME = "C2A-0hz-mirror0.txt"
    # NAME = "C2A-0hz-mirror1.txt"
    # NAME = "C2A-207hz-mirror00.txt"
    # NAME = "C2A-207hz-mirror01.txt"
    # NAME = "C2A-207hz-mirror02.txt"
    # NAME = "C2A-207hz-mirror03.txt"
    # NAME = "C2A-103hz-mirror00.txt"
    # NAME = "C2A-103hz-mirror01.txt"
    # TXT_FILE_PATH = file_str_to_path(NAME, dir=os.path.dirname(__file__)+'\Batch 1')

    ## Batch 4
    # SIGNAL_F = 80.125e6*2 #Hz 
    # SAMPLING_F = 1.0e6 #Hz
    # NAME = "C2-000Hz-0mVpp-00001.txt"
    # NAME = "C2-000Hz-0mVpp-00002.txt"
    # NAME = "C2-143Hz-050mVpp-00000.txt"
    # NAME = "C2-143Hz-050mVpp-00001.txt"
    # NAME = "C2-143Hz-050mVpp-00002.txt"
    # NAME = "C2-143Hz-100mVpp-00003.txt"
    # NAME = "C2-143Hz-100mVpp-00004.txt"
    # NAME = "C2-143Hz-100mVpp-00005.txt"
    # NAME = "C2-143Hz-400mVpp-00002.txt"
    # TXT_FILE_PATH = file_str_to_path(NAME, dir=os.path.dirname(__file__)+'\Batch 4')

    ## Batch 6
    # SIGNAL_F = 80.125e6*2 #Hz 
    # SAMPLING_F = 1.0e6 #Hz
    # NAME = "C2-000Hz-000mVpp-05478Correction00000.txt"
    # NAME = "C2-000Hz-000mVpp-05478Correction00001.txt"
    # NAME = "C2-171Hz-140mVpp-05478Correction00000.txt"
    # NAME = "C2-171Hz-140mVpp-05478Correction00001.txt"
    # TXT_FILE_PATH = file_str_to_path(NAME, dir=os.path.dirname(__file__)+'\Batch 6')

    ## Batch 7a
    SIGNAL_F = 80.125e6*2 #Hz 
    SAMPLING_F = 1.0e6 #Hz
    # NAME = "C1_066Hz_160mVpp_00002.txt"
    NAME = "C1_150Hz_200mVpp_00000.txt"
    TXT_FILE_PATH = file_str_to_path(NAME, dir=os.path.dirname(__file__)+'\Batch 7a')
    FITTING = True
    print(f"{TXT_FILE_PATH = }")

    def get_f_and_A_from_filename(txt):
        Hz_pos = txt.find('Hz')
        mVpp_pos = txt.find('mVp')
        f = int(txt[Hz_pos-3:Hz_pos])
        A = int(txt[mVpp_pos-3:mVpp_pos])
        return f, A
    f, mVpp = get_f_and_A_from_filename(NAME)

    readTxt = fr.read_oscilliscope_txt(TXT_FILE_PATH)
    meta, trace = fr.parse_oscilliscope_txt(readTxt)

    # Prime the algorithm
    N, M = freq_ratio(signal=SIGNAL_F, sample=SAMPLING_F)
    print(f"[Debug] {N = }")
    print(f"sampling interval = {1/SAMPLING_F}")
    print(f'{meta["Sample Interval"][0] = }')

    # Run the phase algorithm
    signal = signal_from_trace(np.asarray(trace))
    phases = signal_to_phase(signal, N, 2*pi/N)
    phases = phase_reconstruction(phases, 5)
    # phases = signal_to_phase(signal, N, 2*pi/N, False)[::N]
    # phases = phase_reconstruction(phases, 4.8)

    t_axis = np.arange(start= 0, 
        stop= (int(meta["Record Length"][0])-N) * meta['Sample Interval'][0], step= meta['Sample Interval'][0])

    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.plot(t_axis, phases, color = 'mediumblue')

    ax.set_ylabel(r'$\phi$/rad', useTex = True)
    ax.set_xlabel(r'$t$/s', useTex = True)

    # Get the fitting functions
    if FITTING:
        import EPstandard
        from scipy import optimize
        def fit_fn(x, A, f, phi, m, c):
            return A * np.sin(2*pi*f*x + phi) + m*x+c
        popt, pcov = optimize.curve_fit(fit_fn, t_axis, phases, p0 = [mVpp * 0.05, f, 0, 0, 0])
        fittings = EPstandard.easy_read_popt_pcov(popt, pcov)
        ax.plot(t_axis, fit_fn(t_axis, *popt), color = 'darkgreen', markersize = 3)

        plt.title(f"Mirror Freq = {f} Hz, Piezo Ampl = {mVpp} mVpp, Signal Freq = {SIGNAL_F:6.0F}Hz, \n(A, $f$, $\phi_0$, m, c) = {fittings}", useTex = True )

    plt.show(block = True)
    # D:\NUS\Semester 7\PC4199 Honours Project\20211020 Mirror Reading\Batch 1 Results

if __name__ == '__main__':
    main()