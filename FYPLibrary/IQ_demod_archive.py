import numpy as np
from numpy import pi as pi

# attempts to do array-level operations
# findings: python is too slow to utilise branchless programming
def phase_reconstruction_2_one(phases, phase_advance, get_deltas=False):
    """TODO: Function descriptor."""
    if type(phases) != np.ndarray:
        phases = np.asarray(phases)
    results = np.asarray([phases[0]])
    phase_advance %= 2*pi
    
    mod_min = lambda a,b: a if abs(a) < abs(b) else b
    mod_min = np.vectorize(mod_min)

    case1 = phases[1:] - phases[:-1]
    case2 = case1 + 2*pi
    r = mod_min(case1, case2)
    results = np.append(results, r)

    if get_deltas:
        return np.asarray(r)
    results = np.cumsum(results)
    return results

def phase_reconstruction_2_two(phases, phase_advance, get_deltas=False):
    """TODO: Function descriptor."""
    if type(phases) != np.ndarray:
        phases = np.asarray(phases)
    results = np.asarray([phases[0]])
    phase_advance %= 2*pi

    case1 = phases[1:] - phases[:-1]
    case2 = case1 + 2*pi
    caseA = abs(case1) < abs(case2)
    caseB = abs(case1) > abs(case2)
    r = caseA * case1 + caseB * case2
    results = np.append(results, r)

    if get_deltas:
        return np.asarray(r)
    results = np.cumsum(results)
    return results

# python for-loop version of existing vectorising IQ method; also too slow
def phase_reconstruction_2a(ph, phase_advance, get_deltas=False):
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
    results = np.cumsum(np.append(results, r))
    return results

# gradient method reconstruction (jy suggestion)
from math import copysign
sign = lambda x: 0 if x == 0 else int(copysign(1,x))
# sign = lambda x: int(x > 0) - int(x < 0)

def phase_reconstruction_3(phases, phase_advance, grad_tolerance= 0.2, 
    twoPi_tolerance= [1.9, 0], step_size_tol= [2,5], debug=(False, 0,0),
    get_deltas = False):
    """Considers relative sign changes, such that sufficiently large and 
    consecutive phase increments that add up to approximately 2pi will 
    be deemed as an errenous 2pi jump to be subtracted away.
    
    Inputs
    ----
    phases (iterable): phase generated from signal_to_phase functions; 
        optional parameter in sig_to_phase **must** be True
    phase-advance (float): subtracts out advancement of wave from 0. 
        usually alr set to 0 as accounted for in `signal_to_phase`
    grad_tolerance (float): phase changes < abs(grad_tol) are considered
        as phase increments too small to be considered in resulting a 
        2pi jump
    twoPi_tolerance (2-elem array): consecutive sign changes that sum to
        absolute value of between 2pi - list[0] and 2pi + list[1]
    step_size_tol (2-elem array) 2 pi jump is considered valid iff 
        number of consecutive steps of sufficiently large size are 
        within this parameter (inclusive)
    debug (3-elem array): verbose log for this function between [1] and
        [2] index specified.
    
    Outputs
    ----
    phase (np iterable): corrected (hopefully) phase"""
    # Variables passed are valid
    # TODO
    
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

# larger window for time average (ari suggestion) (incomplete)
from numpy import average as average
def phase_reconstruction_4(phases, phase_advance, s=100, q=8, debug=(False, 0,0),
    get_deltas = False):
    """Two pi jump is detected by comparing against q steps ago in a 
    buffer of size s. If there is a 2pi jump, the 2pi difference would 
    be observed from q-abit steps, which is almost certainly indicates a 
    definite 2pi jump to be subtracted away.
    
    Inputs
    ----
    phases (iterable): phase generated from signal_to_phase functions; 
        optional parameter in sig_to_phase **must** be True
    phase-advance (float): subtracts out advancement of wave from 0. 
        usually alr set to 0 as accounted for in `signal_to_phase`
    s (int): buffer size
    q (int): number of datapoints back to subtract the phase from to 
        determine if a 2pi jump has occured
    debug (3-elem array): verbose log for this function between [1] and
        [2] index specified.
    
    Outputs
    ----
    phase (np iterable): corrected (hopefully) phase"""
    # Variables passed are valid
    # TODO
    
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
        # TODO
        # if len(hist_timeavg) > s:
        #     hist_timeavg.pop(0)
                
        # if abs(timeavgarr[i+N]-timeavgarr[i]) \approx 2pi for a few i:
        #     increase by 2pi

    # if get_deltas:
    #     return np.asarray(results)
    # results = np.cumsum(results)
    return # np.asarray(hist_timeavg)

# original phase_reconstruction_4
def rmsc(arr):
    # returns rms(arr - arr.mean)
    arr = np.asarray(arr)
    t = arr - arr.mean()
    return np.sqrt(np.mean(np.square(t)))

class fix2pi():
    def __init__(self):
        self.stream = []
        self.shift = 0
        self.shift_counter = 0
        self.prevjump = 0
        self.buflen = 0
        self.cutoff = 64
        self.rollen = 8
        self.rolavg = []
        self.counter = 0
        
        self.tol = 1 # 2 * rms(noise) + a bit more
        
    def process(self,data):
        # main loop that takes in 
        data += self.shift*2*pi
        self.stream.append(data)
        self.counter +=1

        # when the buffer is not filled, populate and terminate fn call
        if self.buflen != self.cutoff:
            self.buflen +=1
            if self.buflen == self.rollen:
                self.rolavg.append(sum(self.stream)/self.rollen)
            elif self.buflen >= self.rollen:
                self.rolavg.append((self.stream[-1]-self.stream[-1-self.rollen])/self.rollen+self.rolavg[-1])
            return None
        # addition of 2 numbers to get rolling average
        self.rolavg.append((self.stream[-1]-self.stream[-1-self.rollen])/self.rollen+self.rolavg[-1])
        
        # check if jump is centered around 2 pi
        jump = self.rolavg[-1] - self.rolavg.pop(0)
        jump_sub = abs(jump) - 2*np.pi
        if abs(jump_sub) < self.tol: # make sure the signs are the same
            if sign(jump) == self.prevjump:
                self.shift_counter += 1
            else:
                self.shift_counter = max(0, self.shift_counter-1)
        else:
            self.shift_counter = max(0, self.shift_counter-1)
        self.prevjump = sign(jump)

        if self.shift_counter == 0: # condition need to change
            # if 0.5*self.tol < rmsc(self.stream) < 1.5*self.tol:
            #     self.tol = 0.8*self.tol + 0.2*rmsc(self.stream)
            # print(f"{self.counter = :>10}, {self.tol = :>7.4f}")
            if self.tol > np.pi:
                print("abvc")
                self.tol = 0
                self.counter -= self.rollen

        # print(f"{jump_sub = }, {self.prevjump = }")
        if self.shift_counter >= self.cutoff - self.rollen -23:
            print(f"jump detected at {self.counter:>10}")
            self.shift_counter = 0
            # implement shift
            self.shift -= self.prevjump
            # change rolling average
            jump_arg = np.argmax([abs(self.stream[i]-self.stream[i-1]) for i in range(1, self.buflen)]) + 1
            for i in range(jump_arg, self.buflen+1): self.stream[i] -= self.prevjump * 2 * pi
            for i in range(max(0,jump_arg-self.rollen),min(jump_arg,self.cutoff-self.rollen+1)): self.rolavg[i] -= self.prevjump*(i-jump_arg)*2*pi/self.rollen
            for i in range(jump_arg,self.cutoff-self.rollen+1): self.rolavg[i] -= self.prevjump*2*pi
            # print(self.prevjump)
        return self.stream.pop(0)
        # self.stream.pop(0)
        # return jump

    def get_buffer(self):
        return self.stream