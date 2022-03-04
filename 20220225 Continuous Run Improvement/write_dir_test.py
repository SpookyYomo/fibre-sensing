# . Honours Module Folder
# ├ FYPLibrary
# | ├ file_reading.py
# | └ IQ_demod.py
# ├ 20210914 Initial Ball Drop
# | └ C2-ball-drop-00002.txt
# └ Sub Project folder
#   └ this script(.py)

# Initialiastion: Directory appending for my system. Vary the directories as necessary.
import sys, os.path

# Add the FYPLibrary folder in the same level as this folder.
if os.path.dirname(os.path.dirname(__file__)) not in sys.path:
    sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'FYPLibrary'))

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

def main():
    files = get_files()
    print(f"Files selected: {files}")
    DIR_WRITE = write_directory(files)
    print(f"Results will be written to: {DIR_WRITE}")
    return

if __name__ == '__main__':
    main()