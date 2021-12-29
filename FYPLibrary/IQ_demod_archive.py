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