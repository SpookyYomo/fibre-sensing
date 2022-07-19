from numpy import gcd

class phasemeter:
    """Fully describe properties of phasemeter experiment."""
    def __init__(self, fcarrier=None, fsample=None,
            N=None, M=None) -> None:
        self.__fcarrier = fcarrier
        self.__fsample = fsample
        self.__N = N
        self.__M = M

    def freq_ratio(carrier, sample):
        """"Determines the simplest ratio M/N between carrier and sample
        frequency.
        
        f_sampler = (N/M) * f_carrier 

        Input
        ----
        carrier (float): frequency of carrier (centre)
        sample (float): frequency of sampling

        Output
        ----
        N (int): Number of sampling points for 2pi advancement in carrier
        M (int): Number of wavelengths of carrier that occur in N periods of 
            sampler
        """
        # Assumes sampling frequency < carrier frequency
        if sample < carrier:
            raise ValueError("Sampling Frequency < Carrier Frequency assumption has been violated!")

        # original use case of sampling frequency being M/N multiple of the 
        # carrier frequency. need to modify accordingly due to turning the 
        # sampling frequency down.
        # consider the time for 1 period of sampler, how many wavelengths of 
        # carrier has occured
        # then determine how many times (N) to then return back to 2pi for 
        # carrier
        # eg: 10+2/3 wavelengths of 160MHz carrier occurs for 15MHz sampler 
        # => N = 3, for which 32 wavelengths of carrier occurs
        carrier, sample = int(carrier), int(sample)
        greatest_common_divisor = gcd(carrier, sample)
        carrier /= greatest_common_divisor
        sample /= greatest_common_divisor
        N, M = int(sample), int(carrier)
        return N, M

    @property
    def N(self):
        return self.__N

    @property
    def M(self):
        return self.M

    @property
    def fCarrier(self):
        return self.__fcarrier

    @property
    def fSample(self):
        return self.__fsample

    @N.setter
    def N(self, Nval):
        if isinstance(self.fCarrier, None) and isinstance(self.fSample, None):
            pass
        else: # Incomplete logic
            pass