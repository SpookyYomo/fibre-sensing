# This script will obtain the laser linewidth.

# . Honours Module Folder
# ├ FYPLibrary
# | ├ file_reading.py
# | └ IQ_demod.py
# └ Experiments 1 (Optional Depth)
#   └ Sub Project folder
#     ├ traces folder
#     | └ laserwidth-1km.dat
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
        "traces", "laserwidth-1km.dat")

    # Step 2: Read appropriate data into table
    data_table = pd.read_csv(file, sep= '\s+', header= None, 
        names= ('freq_v', 'dBm_v'))

    # Step 3: Generate relevant data
    freq = data_table['freq_v'].values
    freq_MHz = freq / 1e6
    dBm = data_table['dBm_v'].values
    yvals = np.power(10, dBm)

    # Step 4: Plot
    matplotlib.rcParams.update({'font.size': 11})
    fig, ax = plt.subplots(nrows=1, ncols=1)
    mu = 141.214
    my_p02 = [245115e-6, 0, -35.9321]
    ax.plot(freq_MHz-mu, my_cauchy2(freq_MHz-mu, *my_p02), color = 'red')
    ax.scatter(freq_MHz-mu, dBm, s= 3.1)
    # ax.set_yscale('log')
    
    ax.set_ylabel(r'Power/dBm', usetex= True)
    ax.set_xlabel(r'$f$ Offset/MHz', usetex= True)
    # geometry package boundary -0.875in*2 horizontally
    # figures set to 0.85 * textwidth
    fig.set_size_inches(0.85*(8.25-0.875*2), 2.5) 
    fig.tight_layout()

    # # Saving figure block
    try:
        wd = write_directory(__file__)
        print('Saving plot...')
        NAME = os.path.basename(file)
        fig.savefig(os.path.join(wd, 'laser_linewidth.png'), 
            format= 'png', dpi= 300, pad_inches = 0)
        print('Saved')
    except FileNotFoundError:
        print(f"[Error] Folder should have been created in initialisation block")
        exit()
    plt.close()

if __name__ == '__main__':
    main()