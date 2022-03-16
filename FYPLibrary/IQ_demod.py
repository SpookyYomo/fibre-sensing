"""IQ-demodulation functions.

This file aims to determine the relative phase change with respect to t=0
for a well specified signal frequency and sampling frequency of ACCURATE
timing, where the phase change in the signal can occur spontaneously.
"""

import sys
import os.path
# from pprint import pprint
# from matplotlib.pyplot import hist

from numpy.lib.function_base import average
if os.path.dirname(os.path.abspath('')) not in sys.path:
    sys.path.append(os.path.dirname(os.path.abspath('')))
# pprint(sys.path)
import file_reading as fr
import numpy as np
from numpy import pi as pi

def freq_ratio(signal, sample):
    """"Determines the simplest ratio M/N between signal and sample freq.
    
    f_sampler = (N/M) * f_signal 

    Input
    ----
    signal (float): frequency of signal (centroid)
    sample (float): frequency of sampling

    Output
    ----
    N (int): Number of sampling points for 2pi advancement in signal
    M (int): Number of wavelengths of signal that occur in N period of 
        sampler
    """
    # Assumes sampling frequency < signal frequency

    # original use case of sampling frequency being M/N multiple of the 
    # signal frequency. need to modify accordingly due to turning the 
    # sampling frequency down.
    # consider the time for 1 period of sampler, how many wavelengths of 
    # signal has occured
    # then determine how many times (N) to then return back to 2pi for 
    # signal
    # eg: 10+2/3 wavelengths of 160MHz signal occurs for 15MHz sampler 
    # => N = 3, for which 32 wavelengths of signal occurs
    signal, sample = int(signal), int(sample)
    greatest_common_divisor = np.gcd(signal, sample)
    signal /= greatest_common_divisor
    sample /= greatest_common_divisor
    N, M = int(sample), int(signal)
    return N, M

def phase_advance(f_signal, f_sample):
    """Obtain DeltaPhi for discrete time interval MOD 2pi.

    For a signal frequency f, and sampling frequency f_s, where
    f_s = (N/M) * f, s.t. M,N integer,
    the phase \Delta\phi by which the signal frequency will advance is 
    given by 
    \Delta\Phi = (2 pi f) * T_s <- time between discrete/sampled steps
    = 2 pi T_s / T = 2 pi M / N *MOD 2PI*

    Input
    ----
    f_signal (float): frequency of signal in hertz
    f_sample (float): frequency of sampling in hertz

    Output
    ----
    dphi (float): advancement in phase MOD 2pi
    """
    return 2 * pi * ((f_signal / f_sample) - (f_signal // f_sample))

# yielding phase
def IQ_signal_to_phase(signal):
    """For N=4, the standard IQ demod is useable here. 
    Obtains only 1 out of every 4 data points"""
    # this aims to check if the generalised case function with N set to 
    # 4 does work as expected or not
    num_elements = len(signal)
    cols = num_elements // 4
    signal = signal[:cols*4]
    signal = signal.reshape((cols, 4))
    Is, Qs = signal[:,0], signal[:,1]
    return np.arctan2(Qs, Is)

def _get_Is_and_Qs(signal, N, dphi):
    """Internal function. For specified signal, N and dphi, return all 
    Is and Qs.
    
    Input
    ----
    signal (np array): 1D array of y_i over time.
    N (int): number of samples involved in IQ valuation
    dphi (float): expected value of advancement iln phase between time 
        intervals. obtained from phase_advance function
    
    Output
    ----
    Is, Qs (np array): I and Q value for each data point, with reference
        to N future data points
    """
    len_signal = len(signal)
    dphi = dphi % (2*pi)

    # lookup tables
    sines = np.sin([i * dphi for i in range(N)])
    cosines = np.cos([i * dphi for i in range(N)])

    # summation process # rewrite for direct multiplication
    # print(f"[Debug] {signal[0:N] = }")
    # Is = (2/N) * np.fromiter([np.dot(signal[j:j+N], sines) for j in range(len_signal-N)], 
    #     dtype= np.float64, count= len_signal-N)
    # Qs = (2/N) * np.fromiter([np.dot(signal[j:j+N], cosines) for j in range(len_signal-N)], 
    #     dtype= np.float64, count= len_signal-N)
    Is = (2/N) * np.convolve(signal, sines[::-1], mode= 'valid')
    Qs = (2/N) * np.convolve(signal, cosines[::-1], mode= 'valid')
    return Is, Qs

def signal_to_phase(signal, N, dphi, phase_advancement_correction):
    """Generalized IQ method to determine the phase from the signal.
    
    Over N consequetive data samples, determine the I and Q values for 
    this interval, to then determine the phase at the start of the start 
    of this interval. Utilises Section 3.1 and 3.2.

    Edit (1 Jan): phase_advancement_corr has is now to be set explicit.
    
    Input
    ----
    signal (np array): 1D array of y_i over time.
    N (int): number of samples involved in IQ valuation
    dphi (float): expected value of advancement in phase between time 
        intervals. obtained from phase_advance function
    phase_advancement_correction (bool): to account for increase in phase
        within a period due to sampling every data point for phase but 
        having a constant IQ axis instead of an IQ axis that advances 
        with the signal to account for appropriate delta

    Output
    ----
    phase (np array): phase of signal starting from Nth data point, with
        value of phase mod 2pi
    """
    Is, Qs = _get_Is_and_Qs(signal, N, dphi)
    result = np.arctan2(Qs, Is)

    # accounting for phase advancement that occurs naturally
    # subtracting the accumulated phase by mod 2pi
    if phase_advancement_correction and not np.isclose(dphi, 2*pi):
        len_signal = len(signal)
        # assumes that N * dphi leads back to 2pi
        subtraction = [j * dphi for j in range(N)] * (len_signal//N + 1)
        subtraction = np.asarray(subtraction[:len_signal-N])
        result = result - subtraction
    
    return result

# debuggering/ML code
def get_R_signal(signal, N, dphi):
    """Determine the R value of the phasor over the entire trace signal.
    
    Input
    ----
    signal (np array): 1D array of y_i over time.
    N (int): number of samples involved in IQ valuation
    dphi (float): expected value of advancement in phase between time 
        intervals. obtained from phase_advance function
    
    Output
    ----
    Rs (np array): sqrt(I^2+Q^2) value for each data point, with 
        reference to N future data points
    """
    Is, Qs = _get_Is_and_Qs(signal, N, dphi)
    result = np.sqrt(np.square(Is) + np.square(Qs))
    return result

def get_Is_Qs(signal, N, dphi):
    """Gets coordiantes in IQ plane.
    
    Input
    ----
    signal (np array): 1D array of y_i over time.
    N (int): number of samples involved in IQ valuation
    dphi (float): expected value of advancement in phase between time 
        intervals. obtained from phase_advance function
    
    Output
    ----
    Is, Qs (np array)
    """
    return _get_Is_and_Qs(signal, N, dphi)

# phase plotting
def phase_reconstruction_naive(ph, tol):
    """NAIVE phase reconstruction after arctan.
    
    If the phase changes too much (ie delta more than tol) will lead to 
    a wrap around.
    
    Inputs
    ----
    ph (np iterable): phase generated from signal_to_phase functions
    tol(float): if change in phase exceeds this value, perform a wrap

    Outputs
    ----
    phase (np iterable): phase wrapped around 2pi accounted phase
    """
    # as taken directly from the MATlab code, determine the 
    # reconstructed phase.
    r = ph.copy()
    # phase reconstruction as done in matlab
    for i in range(1,len(ph)):
        if ph[i] - ph[i-1] > tol:
            r[i:] = r[i:] - np.array([2*pi]*(len(ph)-i))
        elif ph[i] - ph[i-1] < -tol:
            r[i:] = r[i:] + np.array([2*pi]*(len(ph)-i))
    
    return r

# tolerance method
def phase_reconstruction(ph, tol, phase_advance= 0):
    """Phase reconstruction after arctan.
    
    If the phase changes too much (ie delta more than tol) will lead to 
    a wrap around. This function attempts to properly unfold.
    
    Inputs
    ----
    ph (iterable): phase generated from signal_to_phase functions; 
        optional parameter in sig_to_phase **must** be True
    tol (float): if change in phase exceeds this value, perform a wrap
    phase-advance (float): subtracts out advancement of wave from 0. 
        usually alr set to 0 as accounted for in `signal_to_phase`

    Outputs
    ----
    phase (np iterable): phase wrapped around 2pi accounted phase
    """
    if tol > 2*pi: # this function is redundant when tolerance cannot be achieved
        return np.asarray(ph) if not isinstance(ph, np.ndarray) else ph

    r = [ph[0]]
    mod_twopi = 0
    phase_advance %= 2*pi
    for i in range(1, len(ph)):
        if ph[i] - ph[i-1] > tol:
            mod_twopi -= 1
        elif ph[i] - ph[i-1] < -tol:
            mod_twopi += 1
        r.append(ph[i] + 2 * pi * mod_twopi - phase_advance)
    r = np.asarray(r) # faster to append a Py list, convert to np
    return r

# delta tracking method
def phase_reconstruction_2(ph, phase_advance, get_deltas=False):
    """Phase reconstruction after arctan.
    
    Cumulatively sum the deviation of phase from expected phase.

    Inputs
    ----
    ph (iterable): phase generated from signal_to_phase functions; 
        optional parameter in sig_to_phase **must** be True
    tol (float): if change in phase exceeds this value, perform a wrap
    phase-advance (float): subtracts out advancement of wave from 0. 
        usually alr set to 0 as accounted for in `signal_to_phase`

    Outputs
    ----
    phase (np iterable): phase wrapped around 2pi accounted phase"""
    results = np.asarray([ph[0]])
    phase_advance %= 2*pi
    x = np.convolve(ph, [1, -1], mode= 'valid') - phase_advance
    r = 2*pi*(x<-pi) + x
    if get_deltas:
        return r
    results = np.cumsum(np.append(results, r))
    return results