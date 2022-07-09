import sys
import os.path
from pprint import pprint
if os.path.dirname(os.path.abspath('')) not in sys.path:
    sys.path.append(os.path.dirname(os.path.abspath('')))
# pprint(sys.path)
import file_reading as fr
import matplotlib.pyplot as plt
from pathlib import Path

if __name__ == '__main__':
    DIRECTORY = os.path.dirname( __file__ )
    txtfilename = 'truncatedTest.txt'
    TXT_FILE_PATH = os.path.join(DIRECTORY, txtfilename)
    TXT_FILE_PATH = Path("D:\\NUS\\Semester 7\\PC4199 Honours Project\\20210914 Oscilliscope Reading Attempt 2\\Nis5data\\C280.020Mhz100kS-nothing00001.txt")
    readTxt = fr.read_oscilliscope_txt(TXT_FILE_PATH)
    print('Read')
    # pprint(readTxt)
    
    meta, trace = fr.parse_oscilliscope_txt(readTxt)
    # print('\nParse')
    # pprint(meta)
    # pprint(trace)

    # trace = fr.trace_error_correct(meta, trace)
    # print('\nCorrection')
    # pprint(trace)

    signal = fr.signal_from_trace(trace)
    plt.plot(signal)
    plt.show()