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
    # print(f"[Debug] {N = }")
    # print(f"[Debug] {len(sines) = }")

    # summation process # rewrite for direct multiplication
    # Is = (2/N) * np.fromiter([np.sum([signal[j+i]*sines[i] for i in range(N)]) for j in range(len_signal-N)], dtype=np.float64)
    # Qs = (2/N) * np.fromiter([np.sum([signal[j+i]*cosines[i] for i in range(N)]) for j in range(len_signal-N)], dtype=np.float64)
    # print(f"[Debug] {signal[0:N] = }")
    Is = (2/N) * np.fromiter([np.dot(signal[j:j+N], sines) for j in range(len_signal-N)], 
        dtype=np.float64)
    Qs = (2/N) * np.fromiter([np.dot(signal[j:j+N], cosines) for j in range(len_signal-N)], 
        dtype=np.float64)
    return Is, Qs

def signal_to_phase(signal, N, dphi, phase_advancement_correction= False):
    """Generalized IQ method to determine the phase from the signal.
    
    Over N consequetive data samples, determine the I and Q values for 
    this interval, to then determine the phase at the start of the start 
    of this interval. Utilises Section 3.1 and 3.2.

    Edit (27 Dec): phase_advancement_corr has been set to false with 
    increasing preference for phase accumulation method in phase 
    reconstruction. If any old graphs have broken, please set the 
    opstional parameter to true (particularly if using 
    phase_reconstruction_1).
    
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

### Experimental
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
    # for i in range(1, len(ph)):
    #     expected_phase = ph[i-1] + phase_advance
    #     # choose between the least difference (especially when expected 
    #     # phase > 2pi)
    #     delta = case1 if abs(case1:= ph[i]-expected_phase) < \
    #         abs(case2:= case1+2*pi) else case2
    #     results.append(delta)
    r = np.fromiter([case1 if abs(case1:= ph[i]-ph[i-1]-phase_advance) < \
            abs(case2:= case1+2*pi) else case2 for i in range(1,len(ph))], \
                dtype= np.float32, count= len(ph)-1)

    if get_deltas:
        return r
    results = np.append(results, r)
    return results


# gradient method reconstruction (jy suggestion)
from math import copysign
sign = lambda x: 0 if x == 0 else int(copysign(1,x))
# sign = lambda x: int(x > 0) - int(x < 0)

def phase_reconstruction_3(phases, phase_advance, grad_tolerance= 0.2, 
    twoPi_tolerance= [1.9, 0], step_size_tol= [2,5], debug=(False, 0,0),
    get_deltas = False):
    """TODO: Function descriptor."""
    # Variables passed are valid
    
    # Initialisation
    results =  [phases[0]]
    phase_advance %= 2*pi
    m, prev_m, css, prev_css = 0, 0, 0, 0 # gradient direction (with diffused 0); number of consecutive same sign
    
    for i in range(1, len(phases)):
        expected_phase = phases[i-1] + phase_advance
        delta = case1 if abs(case1:= phases[i]-expected_phase) < abs(case2:=case1+2*pi) else case2
        prev_m = m
        if abs(delta) > grad_tolerance:
            # "add" 1 if same sign 
            m = m + sign(delta) if sign(m)*sign(delta) == 1 else sign(delta)
        else:
            m = 0
        prev_css = css
        css = css+sign(delta) if sign(delta)==sign(css) else sign(delta)
        
        if debug[0] and debug[2] - debug[1] < 101 and debug[1] <= i <= debug[2]:
            print(f"\n{i = }\n{delta = }")
            print(f"{m = }, {prev_m = },{sign(m)!=sign(prev_m) = }")
            print(f"{step_size_tol[0] <= abs(prev_m) <= step_size_tol[1] = }")
            print(f"{css = }, {sum(results[-prev_css:]) = }")
            print(f"{(2*pi - twoPi_tolerance[0]) <= abs(sum(results[-prev_css:])) < (2*pi + twoPi_tolerance[1]) = }")

        # if m==0 and abs(prev_m) in step_size_tol or sum(results[-prev_css:]) is close to 2*pi:
        if sign(m)!=sign(prev_m) and step_size_tol[0] <= abs(prev_m) <= step_size_tol[1] and \
            (2*pi - twoPi_tolerance[0]) <= abs(sum(results[-prev_css:])) < (2*pi+twoPi_tolerance[1]):
            print(f"Step {i = } has injected {-sign(prev_m)} * 2*pi.")
            for j in range(s:= step_size_tol[0]):
                results[-j] += -sign(prev_m) *2*pi/s
            
        results.append(delta)
    if get_deltas:
        return np.asarray(results)
    results = np.cumsum(results)
    return results

# larger window for time average (ari suggestion)
def phase_reconstruction_4(phases, phase_advance, s=100, q=8, debug=(False, 0,0),
    get_deltas = False):
    """TODO: Function descriptor."""
    # Variables passed are valid
    
    # Initialisation
    results =  [phases[0]]
    phase_advance %= 2*pi
    
    for i in range(1, s):
        expected_phase = phases[i-1] + phase_advance
        delta = case1 if abs(case1:= phases[i]-expected_phase) < abs(case2:=case1+2*pi) else case2

        results.append(delta)

    hist_timeavg = []
    for i in range(s, len(phases)):
        expected_phase = phases[i-1] + phase_advance
        delta = case1 if abs(case1:= phases[i]-expected_phase) < abs(case2:=case1+2*pi) else case2

        results.append(delta)

        hist_timeavg.append(timeavg_recent := average(phases[-s:]))
        # if len(hist_timeavg) > s:
        #     hist_timeavg.pop(0)
                
        # if abs(timeavgarr[i+N]-timeavgarr[i]) \approx 2pi for a few i:
        #     increase by 2pi

    # if get_deltas:
    #     return np.asarray(results)
    # results = np.cumsum(results)
    return np.asarray(hist_timeavg)