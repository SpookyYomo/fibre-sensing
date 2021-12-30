# Initialiastion: Directory appending for my system. Vary the directories as necessary.
import sys
import os.path
# Add the FYPLibrary folder in the same level as this folder.
if os.path.dirname(os.path.dirname(__file__)) not in sys.path:
    # print(f"{os.path.dirname(os.path.dirname(__file__)) = }")
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'FYPLibrary'))

from pprint import pprint
# pprint(sys.path)

# Import Modules
from file_reading import *

def main():
    # NAME = "C1_150Hz_200mVpp_00000.txt"
    TXT_FILE_PATH = r'D:\NUS\Semester 7\PC4199 Honours Project\20211020 Mirror Reading\Batch 7a\C1_066Hz_160mVpp_00002.txt'
    # TXT_FILE_PATH = file_str_to_path(NAME, dir=os.path.dirname(__file__)+'\Batch 7a')
    meta, trace = parse_and_read_oscilliscope_txt(TXT_FILE_PATH)
    pprint(meta)
    print(f"{meta['Record Length'] = }")

if __name__ == '__main__':
    main()