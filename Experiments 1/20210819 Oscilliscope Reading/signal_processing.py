import numpy as np

def discrete_fourier(signal):
    """
    Obtains the fourier transform of a discretized EVENLY spaced signal.

    Input
    -----
    Signal: Array_like of (real/complex)

    Output
    -----
    FFT(Signal) 
    """

    return np.fft.fft(signal)

if __name__ == '__main__':
    import file_reading as fr
    from pathlib import Path
    readTxt = fr.read_oscilliscope_txt(Path('D:\\NUS\\Semester 7\\PC4199 Honours Project\\20210812 Oscilliscope Matters\\Test Files\\C2Test16Aug202100001.txt'))
    meta, trace = fr.parse_oscilliscope_txt(readTxt)
    trace = fr.trace_error_correct(meta, trace)
    signal = fr.signal_from_trace(trace)
    del readTxt

    import matplotlib.pyplot as plt
    freq = discrete_fourier(signal)
    plt.plot(freq)
    plt.show()