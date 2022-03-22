# This script will obtain the AOM diffraction effeciency graph.

# . Honours Module Folder
# ├ FYPLibrary
# | ├ file_reading.py
# | └ IQ_demod.py
# └ Sub Project folder
#   ├ traces folder
#   | └ AOM-diffraction-eff.dat
#   └ this script(.py)

# Initialiastion: Directory appending for my system. Vary the directories as necessary.
import sys, os.path

# Add the FYPLibrary folder in the same level as this folder.
if os.path.dirname(os.path.dirname(__file__)) not in sys.path:
    sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'FYPLibrary'))

# from pprint import pprint
# pprint(sys.path)

# within this jupyter notebook
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

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

def my_cauchy2(x, gamma, a , b):
    return b+10*np.log10(gamma**2/((x-a)**2+gamma**2))

def main():
    # Step 1: Get file
    file = os.path.join(os.path.dirname(__file__),
        "traces", "AOM-diffraction-eff.dat")

    # Step 2: Read appropriate data into table
    data_table = pd.read_csv(file, sep='\s+')

    # Step 3: Generate relevant data
    freq = data_table['aom-freq'].values
    eff = data_table['1st_order_powermeter'].values

    # Step 4: Plot
    matplotlib.rcParams.update({'font.size': 11})
    fig, ax = plt.subplots(nrows=1, ncols=1)
    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.plot(freq, eff, color = 'mediumblue', linewidth= 1)
    ax.set_ylabel(r'Diffraction Efficiency/\%', usetex= True)
    ax.set_xlabel(r'$f$/MHz', usetex= True)
    ax.set_ylim(20,80)

    # geometry package boundary -0.875in*2 horizontally
    # figures set to 0.85 * textwidth
    fig.set_size_inches(0.85*(8.25-0.875*2), 2.5) 
    fig.tight_layout()

    # # Saving figure block
    try:
        wd = write_directory(__file__)
        print('Saving plot...')
        NAME = os.path.basename(file)
        fig.savefig(os.path.join(wd, 'diffraction-efficiency.png'), 
            format= 'png', dpi= 300, pad_inches= 0)
        print('Saved')
    except FileNotFoundError:
        print(f"[Error] Folder should have been created in initialisation block")
        exit()
    plt.close()

if __name__ == '__main__':
    main()