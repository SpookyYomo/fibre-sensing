"""File reading functions.

This file aims to help with the reading of LeCroy output files.
"""

import os, re
from signal import signal

# from numpy.core.defchararray import array
import numpy as np
import pandas as pd
import lecroyparser as lcp

### General ###
def transpose(m):
    """Transposes a matrix."""
    return list(map(list,zip(*m)))

def table_trim(m, r):
    """Cuts off the title columns and blank space, and returns the numpy-ed array."""
    # Yield first empty row
    empty = ""
    if r == None:
        try:
            index_first_empty = m[0].index(empty)
        except:
            index_first_empty = 0
            print('Data table beginning from row 0!')
    else:
        index_first_empty = r - 2   

    # Cast to numpy array to return
    result = []
    for i in range(len(m)):
        hold = m[i][index_first_empty+1:]
        result.append(np.array(hold))
    return result

def file_str_to_path(filename, dir = os.path.dirname(__file__)):
    """From file name, get a path object.
    
    Inputs
    ------
    filename: String of file, with .txt;
    dir: os.path object that shows folder of file. Optional.
    """
    return os.path.join(dir, filename)

### Data Reading functions ###
def line_splitter(line):
    """Splits a line into array of elements without splitting substrings
    Inputs: string
    Outputs: array

    Alternatively, this can be done by shlex.split module, but is not 
    used due to backward compatibility requirement.
    STRIPS OUT THE DOUBLE QUOTES
    """

    if '"' not in line:
        result = line.split()
    else:
        # result = [p for p in re.split("( |\\\".*?\\\"|'.*?')", line) if p.strip()]
        # result = [p if '"' not in p else p.strip('"') for p in result]
        result = [p if '"' not in p else p.strip('"') for p in re.split("( |\\\".*?\\\"|'.*?')", line) if p.strip()]
        # result = [p for p in re.split("( |\\\".*?\\\"|'.*?')", line) if p.strip('"')]
    return result

def read_oscilliscope_txt(txtfilepath):
    """DEPRECIATED Use: parse_and_read_oscilliscope_txt(...).
    From text file, output into python array.
    
    Inputs
    ------
    txtfilepath: os.path object

    Outputs
    -------
    python array of floated (if possible) data
    """
    rows = []

    # `with` block ensures file is closed when complete
    with open(txtfilepath, 'r') as txtfile:
        # file is large, attempt to obtain line by line
        txtlines = txtfile.readlines()
        for line in txtlines:
            # removes spacebar, gets individual items in the line, remove \n in line
            hold = []
            line = line[:-1]

            # split either by special or default
            if '"' not in line:
                items = line.split() # reduce function calling
            else:
                items = line_splitter(line)
            
            # everything so far is still in string format
            for item in items:
                try:
                    hold.append(float(item))
                except ValueError: # Floats a non-floatable
                    hold.append(item)

            if hold != []:
                rows.append(hold)
    return rows

def parse_oscilliscope_txt(txtarray):
    """DEPRECIATED Use: parse_and_read_oscilliscope_txt(...).
    Returns a proper data structure of the txt array.
    
    Inputs
    ------
    txtarray: Python Array
        Array that has converted cell text lines to individual elements 
        of python array.

    Outputs
    -------
    metadata: Dictionary
        Non-oscilliscope data that was saved in the txt
    trace: Python array (n by 2)
        n rows of 2-element tuple: (timing, oscilliscope trace). 
    """

    # with reference to the txt structure, it appears that the original 
    # format has 5 columns of data, with the first 3 going to metadata
    # given the nature of the text format, and spacebars having been 
    # reomoved when converting into python array, we shall assume that 
    # if the row is 5 columns wide, it contains meta data in the first 3 
    # elements, and if two, contains only trace data
    metadata = {}
    trace = []

    for line in txtarray:
        if len(line) == 2: # Hard-coded
            # first statement to reduce assembly point jump
            trace.append(tuple(line))
        elif len(line) == 5: # Hard-coded
            metadata[line[0]] = tuple(line[1:3]) # key = value
            trace.append(tuple(line[3:]))
        else:
            raise ValueError('Unrecognised data structure!')
            print(f"Error encountered on line [{txtarray.index(line)+1}]: {line}")

    # trace = np.asarray(trace) # floating conversion is buggy

    return metadata, trace

def parse_and_read_oscilliscope_txt(txtfilepath):
    """From text file (ASCII w/ header), output into meta-data and trace.
    
    Inputs
    ------
    txtfilepath: os.path object
    
    Outputs
    -------
    metadata: Dictionary
        Non-oscilliscope data that was saved in the txt
    trace: (n by 2) numpy array
        n rows of 2-element tuple: (timing, oscilliscope trace). 
    """
    if txtfilepath[-4:] != '.txt':
        raise ValueError(f'Incorrect file type has been used (txt expected): {txtfilepath}')
    # parse by pandas
    trace = pd.read_csv(txtfilepath, delimiter=' ', header=None, 
        usecols=[3,4])
    print('[Notice] Ignore the following mixed types warning.')
    meta_temp = pd.read_csv(txtfilepath, delimiter=' ', header=None,
        usecols=[0,1,2], dtype={'a':str, 'b': np.float32, 'c': str},
        keep_default_na=True).dropna()

    # convert pandas to numpy array
    trace = trace.to_numpy()
    # clean metadata to dictionary
    keys, numericals, units = meta_temp[0], meta_temp[1], meta_temp[2]
    vals = zip(numericals, units)
    meta = dict(zip(keys, vals))
    print(f"[Log] file supplied = {txtfilepath}\n[Log] {meta = }")
    return meta, trace

def parse_and_read_oscilliscope_trc(trcfilepath):
    """From Le Croy trc file (binary w/ word), output into numpy-trace.
    
    Inputs
    ------
    txtfilepath: os.path object
    
    Outputs
    -------
    time-axis: (n by 1) numpy array of timing
    signal: (n by 1) numpy array of signal
    """
    if trcfilepath[-4:] != '.trc':
        raise ValueError(f'Incorrect file type has been used (trc expected): {trcfilepath}')
    
    scope_data = lcp.ScopeData(path= trcfilepath)
    time, *signal = scope_data.parseFile(path= trcfilepath)
    print(f"[Log] file supplied = {trcfilepath}\n[Log] meta is not embedded")
    return time, signal

def timeaxis_from_trace(trace):
    """Extract temporal component of trace structure.

    Input
    -----
    Trace: (n x 2) python/numpy array

    Output
    -----
    timeaxis: numpy array (n elements)

    """
    # return np.fromiter(map(lambda x: x[1], trace), dtype=np.float64)
    if not isinstance(trace, np.ndarray): #if type(trace) != np.ndarray:
        trace = np.asarray(trace)

    # asserts to be of 2 dimensions with 2 element per row
    assert trace.ndim == 2
    assert trace.shape[1] == 2 

    return trace[:,0]

def signal_from_trace(trace):
    """Extract signal component of trace structure.

    With the metadata and trace corrected, we assumed that regular 
    spacing has since been guranteed. Now obtain just the signal to 
    perform FFT to determine the fourier signals.
    
    Input
    -----
    Trace: (n x 2) python/numpy array

    Output
    -----
    Signal: numpy array (n elements)

    """
    # return np.fromiter(map(lambda x: x[1], trace), dtype=np.float64)
    if not isinstance(trace, np.ndarray): #if type(trace) != np.ndarray:
        trace = np.asarray(trace)

    # asserts to be of 2 dimensions with 2 element per row
    assert trace.ndim == 2
    assert trace.shape[1] == 2 

    return trace[:,1]

### Text parser ###
def _text_trimmer(elem):
    """Removes only alphabets from a string."""
    return ''.join(i for i in elem if not i.isalpha())

def numerical_variables_from_name(name, delimiter = '-'):
    """Removes alphabets to get numbers from Lecroy filename.
    
    Example: 
    'C2-341Hz-250mVpp-040Correction00003.txt' -> [2, 341, 250, 40]
    'C1-M150Hz-M200mVpp--050HzDelta00000.txt' -> [1, 150, 200, -50]
    'C1-M150Hz-M200mVpp-80.125050Hz00000.txt' -> [1, 150, 200, 80.125050]

    Input
    -----
    name: string, as given

    Output
    -----
    result: list of floats of variables extracted
    """
    name = name[:-9] # Remove 00000.txt
    name = name.split(delimiter)
    if name[-1] == '': name.pop() # Special case C2-text-00000.txt

    if delimiter == '-' and '' in name:
        # Reinsert - sign to next element in list since removed by .split
        indices = [i for i, x in enumerate(name) if x=='']
        for i in indices:
            name[i+1] = '-' + name[i+1]
    # Strip text
    name = map(_text_trimmer, name)

    # Remove blank elements
    name = filter(lambda x: x != '', name)
    # for _ in filter(lambda x: x == '', name):
    #     name.remove('')

    # Yield numbers only
    result = list(map(float, name))
    return result

### Debugging

def trace_error_correct(metadata, trace):
    """determine if trace has unsual time stamping, guided by meta-data
    
    notifies if only the last trace data point is wrong, and will correct the time stamp accordingly.
    if more than one timing data stamp is wrong, halts the programme.
    
    generally not used because oscilliscope timing output is unstable

    Inputs
    ------
    metadata: Python dictionary
        supplied out from parse_oscilliscope_txt function.
        minimally contains keys "Record Length", "Sample Interval" and, 
        "Horizontal Offset"
    trace: Iterable of 2 elements
        supplied out from parse_oscilliscope_txt function.
        contains (time, value) per iterable

    Outputs
    -------
    trace: same as input, but with fixed (or nothing changed) trace
    """
    
    ## metadata strucure assumed to contain
    # 1. "Record Length": number of datapoints (and hence nyumber of 
    #   rows) 
    # 2. "Sample Interval": time stamp between one data point to the 
    #   next + units of measurement
    # 3. "Horizontal Offset": first data point (?) + units of measurement

    # to avoid assumption that the first data point necessarily is given 
    # by horizontal offset, we just check instead if the linear spacing 
    # as given is indeed correct, but just notify that the first data 
    # point does not agree

    # notifies that meta data does not align with first data point
    if metadata["Horizontal Offset"][0] != trace[0][0]:
        print("[Warning] Meta-data first data point does not match trace timing.")

    # check that all time stamps are equal up to 2nd last data point
    # breaks if found to have any differing
    if metadata["Record Length"][1] != 'Points':
        print(f"{metadata['Record Length'][1] = }")
        raise ValueError('Unrecognised Metadata structure for "Record Length"')
    if metadata["Horizontal Offset"][1] != metadata["Sample Interval"][1]:
        print(f"{metadata['Horizontal Offset'][1] = }")
        raise ValueError('Different units of time')
    
    n = int(metadata["Record Length"][0]) # number of samples
    t0 = trace[0][0] # intial timing, as obtained from trace
    delta = metadata["Sample Interval"][0] 
    reltol = 0.08
    abstol = delta/100 # hardcoded

    data_time = np.fromiter(map(lambda x: x[0], trace[:-1]), dtype=np.float64)
    differences = data_time[1:-1] - data_time[:-2]
    all_match = np.allclose(100*differences, [100*delta]*(n-3), rtol=reltol, atol= 100*abstol)
    if not all_match:
        print(f"[Notice] Mismatch detected! Determining invalid timing.")
        for j in range(n-3):
            if not np.isclose(differences[j], delta, rtol=reltol, atol= abstol):
                print(f"[Error] First error encountered on row {j+1}. Trace value found to be \n {trace[j-1:j+2] = }\nThis is in contrast to the expected {[t0 + i*delta for i in range(j-1, j+2)] =}")
                raise ValueError('Timing data for trace found to have been incorrect')
        print(f"[Notice] Numpy np.allclose bug found to be faulty.")
    
    
    # data_time = np.fromiter(map(lambda x: x[0], trace[:-1]), dtype=np.float64)
    # expected_time = np.fromiter(map(lambda i: t0 + i*delta, range(n-1)), dtype=np.float64)
    # # if any([trace[i][0] != t0 + i*delta for i in range(n-1)]):
    # if not np.allclose(data_time, expected_time, atol = abstol):
    #     print(f"[Notice] Mismatch detected! Determining invalid timing.")
        
    #     matching = np.fromiter(map(lambda j: 
    #         np.isclose(data_time[j], expected_time[j], atol = abstol), range(n-1)), dtype=bool)
    #     if False in matching:
    #         js = np.where(matching == False)
    #         for j in js:
    #             print(f"[Error] First error encountered on row {j+1}. Trace value found to be \n {trace[j-1:j+2] = }\nThis is in contrast to the expected {[t0 + i*delta for i in range(j-1, j+2)] =}")
    #         raise ValueError('Timing data for trace found to have been incorrect')
    #     else:
    #         print(f"[Notice] Numpy np.allclose bug found to be faulty.")
        # flag, j = False, 0 # initiallisation value to determine row with error
        # while not flag and j < n:
        #     if not np.isclose(trace[j][0], t0 + j*delta, abstol):
        #         flag = True
        #         print(f"[Error] First error encountered on row {j+1}. Trace value found to be \n {trace[j-1:j+2] = }\nThis is in contrast to the expected {[t0 + i*delta for i in range(j-1, j+2)] =}")
        #     j += 1

        # raise ValueError('Timing data for trace found to have been incorrect')
    
    # correct the last data point if necessary
    if trace[-1][0] != t0 + (n-1)*delta:
        print(f"[Notice] Due to oscilliscope behaviour, the last timing trace was found to be wrong.\nThis has been accounted for and will be edited.\n\nOld value: {trace[n-1][0] = }")
        # assumes data structure as given in input of this function
        trace[-1] = (t0 + (n-1)*delta, trace[-1][1])
        print(f"New value: {trace[n-1][0] = }")

    return trace