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

# Imports for this script
from os import listdir
from file_reading import *
from IQ_demod import *
import matplotlib.pyplot as plt
from math import pi

# Yield all files in batch folder
folderpath = os.path.join(os.path.dirname(__file__),'traces')
for i, f in enumerate(files:= listdir(path= folderpath)):
    print(f"{i:>3}: {f}")
# my_input = int(input("Please specify which file you would like to generate a periodogram of: "))
my_input = 1
NAME = files[my_input]
TXT_FILE_PATH = file_str_to_path(NAME, dir=folderpath)
print(f"{TXT_FILE_PATH = }")
_, mirror_f, mVpp, delta_corr = numerical_variables_from_name(NAME)
print(f"{mirror_f = }, {mVpp = }, {delta_corr = }")

zoom_in_start = 0
zoom_in_end = zoom_in_start + 500000

SIGNAL_F = 80.125e6*2 #Hz 
SAMPLING_F = 1.0e6 #Hz
N, _ = freq_ratio(signal=SIGNAL_F, sample=SAMPLING_F)

meta, trace = parse_and_read_oscilliscope_txt(TXT_FILE_PATH)
signal = signal_from_trace(np.asarray(trace))

# phase reconstuction block
phases = signal_to_phase(signal, N, 2*pi/N, phase_advancement_correction= False)
t_axis = np.arange(start= 0, 
    stop= (int(meta["Record Length"][0])-N+1) * meta['Sample Interval'][0],
    step= meta['Sample Interval'][0])
phases = phase_reconstruction_3(phases, 2*pi/N, grad_tolerance=0.3)

fig, ax = plt.subplots(nrows=1, ncols=1)
ax.plot(t_axis, phases[zoom_in_start:zoom_in_end], color = 'mediumblue')
ax.set_ylabel(r'$\phi$/rad', useTex = True)
ax.set_xlabel(r'$t$/s', useTex = True)