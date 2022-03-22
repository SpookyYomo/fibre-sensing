# This script will obtain the std deviation graph from processed logs.

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

# within this jupyter notebook
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.dates as mdates

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
    # Step 1: Get file
    file = os.path.join(os.path.dirname(__file__),
        "traces", "processedLogs-20220224-1557.txt")

    # Step 2: Read appropriate data into table
    data_table = pd.read_csv(file, sep= '\t', header= None, 
        names = ('id', 'mean_v', 'std_dev', 'max_v', 'min_v'))

    # Step 3: Generate x-axis data
    starttime = np.datetime64('2022-02-24 15:57')
    my_time_delta = data_table['id'].values * 61
    time_axis = starttime + pd.to_timedelta(my_time_delta, 'S') 
    y_vals = data_table['std_dev'].values

    # Step 4: Plot
    matplotlib.rcParams.update({'font.size': 11})
    fig, ax = plt.subplots(nrows=1, ncols=1)
    locator = mdates.AutoDateLocator()
    formatter = mdates.ConciseDateFormatter(locator)
    formatter.formats = ['%y',  # ticks are mostly years
                            '%b',       # ticks are mostly months
                            '%d',       # ticks are mostly days
                            '%H:%M',    # hrs
                            '%H:%M',    # min
                            '%S.%f', ]  # secs
    # these are mostly just the level above...
    formatter.zero_formats = [''] + formatter.formats[:-1]
    # ...except for ticks that are mostly hours, then it is nice to have
    # month-day:
    formatter.zero_formats[3] = '%d-%b'

    formatter.offset_formats = ['',
                                '%Y',
                                '%b %Y',
                                '%d %b %Y',
                                '%d %b %Y',
                                '%d %b %Y %H:%M', ]
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)

    ax.bar(time_axis, y_vals, width= 0.0008)
    ax.set_ylabel('Standard Deviation $\sigma$/rad', usetex= True)
    ax.set_xlabel('Time', usetex= False)
    # geometry package boundary -0.875in*2 horizontally
    # figures set to 0.85 * textwidth
    fig.set_size_inches(0.85*(8.25-0.875*2), 2.50) 
    fig.tight_layout()

    # # Saving figure block
    try:
        wd = write_directory(__file__)
        print('Saving plot...')
        NAME = os.path.basename(file)
        fig.savefig(os.path.join(wd, 'std_dev_' + NAME[:-4] + '.png'), 
            format= 'png', dpi= 300, pad_inches = 0.01)
        print('Saved')
    except FileNotFoundError:
        print(f"[Error] Folder should have been created in initialisation block")
        exit()
    plt.close()

if __name__ == '__main__':
    main()