import os
import csv
import numpy as np

### Data processing ###
def fMultiply(x, y, dx, dy):
    """Returns z = x*y, dz. Assumes dx, dy > 0"""
    z = x * y
    dz = np.abs(x) * dy + dx * np.abs(y)
    return z, dz

def fDivision(x, y, dx, dy):
    """Returns z = x/y, dz"""
    z = x/y
    # if dx == 0:
    #     dz = np.abs(z) * np.abs(dy/y)
    # elif dy == 0:
    #     dz = np.abs(z) * np.abs(dx/x)
    # else:
    dz = np.abs(z * np.sqrt( (dx/x)**2.0 + (dy/y)**2.0 ))
    return z, dz

def fLog(x, dx):
    """Returns z = log(x), dz. Assumes dx > 0"""
    z = np.log(x)
    dz = dx / np.abs(x)
    return z, dz

def gen_error(x, dx):
    """For a measured variable `x` with all values having the same error
    `dx`, generate the appropriate dx array."""
    return np.ones_like(x) * dx

### Data Reading functions ###
def read_ExperimentCSV(csvfilepath):
    rows = ()
    
    with open(csvfilepath) as csvfile:
        file_reader = csv.reader(csvfile)
        for row in file_reader:
            hold = []
            for item in row:
                try:
                    hold.append(float(item))
                except:
                    hold.append(item)
            rows += (tuple(hold), )
    return rows

def transpose(m):
    """ Transposes a matrix."""
    return list(map(list,zip(*m)))

def table_trim(m, r):
    """ Cuts off the title columns and blank space, and returns the numpy-ed array."""
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

def csv_file_to_variable(csvfilename, dir = os.path.dirname(__file__), start_row = None):
    """Specify path and first non-empty row of data file to convert to
    numpy array.
    
    Inputs
    ------
    csvfilename: String of csvfile, with .csv;
    dir: os.path object that shows folder of csvfile. Optional;
    start_row: First non-empty row of csvfile, starting from 1.
    """
    data = read_ExperimentCSV(os.path.join(dir, csvfilename))
    return table_trim(transpose(data), start_row)

### Curve fitting (OLS method) ###
# def fitting_equation(x, m, c):
#     # Please edit the parameters for fitting above and the equation below, for scipy.curvefit!!
#     # The equation is given by y= m*x+c
#     return m*x + c

def get_std_dev(pcov):
    # print(f"The Error found to be:{[np.sqrt(i) for i in pcov[j][j] for j in range(len(pcov))]}")
    return [np.sqrt(pcov[j][j]) for j in range(len(pcov))]

# Curve fitting
# popt, pcov = curve_fit(fitting_equation, x_data, y_data, sigma=resolution_error)

def pretty_print_fitted_vals(popt, pcov):
    """Does not return anything, pretty prints popt, pcov from scipy.curvefit"""
    print(f"Fitting values: \n{popt}\nwith errors:\n{get_std_dev(pcov)}\n")

### Curve fitting (ODR method) ###
from scipy.odr import ODR, Model, Data, RealData
from numpy import power

def ODRf(fittype):
    """Generalisation of ODRlinear etc."""

    def ODRlinear(B, x):
        """Linear function y = m*x+c, for ODR.

        Inputs
        ------
        B: n-array, a vector of fitting parameters.
        x: x is an array of current x-values. same format as x passed to Data or
            real data.
        """
        return B[0]*x + B[1]

    def ODRpower2(B, x):
        """Power function y = a * x**b, for ODR.

        Inputs
        ------
        B: n-array, a vector of fitting parameters.
        x: x is an array of current x-values. same format as x passed to Data or
            real data.
        """
        return B[0]*(x**B[1])

    def ODRpower3(B, x):
        """Power function y = a * x**b + c, for ODR.

        Inputs
        ------
        B: n-array, a vector of fitting parameters.
        x: x is an array of current x-values. same format as x passed to Data or
            real data.
        """
        return B[0]*(x**B[1]) + B[2]
    
    def ODRgaussian(B, x):
        """Gaussian function y = a * exp(-(x-b)^2/2c^2), for ODR.

        Inputs
        ------
        B: n-array, a vector of fitting parameters.
        x: x is an array of current x-values. same format as x passed to Data or
            real data.
        """
        return B[0]* np.exp(- (x-B[1])**2. /B[2] )

    if fittype == 'linear':
        return ODRlinear
    elif fittype == 'power2':
        return ODRpower2
    elif fittype == 'power3':
        return ODRpower3
    elif fittype == 'gaussian':
        return ODRgaussian
    else:
        raise ValueError("fittype specified is unrecognised!")

def ODR_popt_pcov(x, y, sx, sy, ODRmodel, initial_estimate):
    """Use scipy.odr to perform EXPLICIT orthogonal distance regression.
    Obtains the equivalent to scipy.curve_fit's popt and pcov.

    Inputs
    ------
    x: x-axis values
    y: y-axis values
    sx: (standard) uncertainty of x-values
    sy: (standard) uncertainty of y-values
    ODRmodel: a function of form f(B,x), where B is fitting parameters to
        be obtained, and x is the x-axis values
    initial_estimate: initial estimate for B
    """

    my_data = Data(x, y, wd=1./power(sx,2), we=1./power(sy,2))
    my_model = Model(ODRmodel)
    my_odr = ODR(my_data, my_model, beta0=initial_estimate)
    my_odr.set_job(fit_type=0) # 2: OLS, 0: Explicit ODR, 1: Implicit ODR
    my_output = my_odr.run()

    popt = my_output.beta
    pcov = my_output.cov_beta

    return popt, pcov

### Easy popt, pcov reading ###
def easy_read_popt_pcov(popt, pcov, verbose= False):
    """ONLY Prints B[i] = 1.29(36) formatting for all i.
    Assumes? std_dev has a smaller order of magntude than popt."""

    # convert pcov to to cov to std_dev 
    perr= np.sqrt(np.diag(pcov))
    if __name__ == '__main__':
        perr= np.sqrt(pcov)

    if verbose:
        print("Full structure")
        for index, (main_num, error_num) in enumerate(zip(popt, perr)):
            print(f"B[{index}]    = {main_num:30.10f}")
            print(f"B[{index}]err = {error_num:30.10f}")
        print()

    def get_first_two_non_zero_digits(f):
        """ get digit and digits placement"""
        f = float(f)
        s = str("%.17f" % f) + '0' # removes scientific notation of floating
        l = list(s)
        # for index, char in enumerate(s):
        #     if char == '.':
        #         dot_index = index
        #         break
        dot_index = s.find('.')
        for index, char in enumerate(s):
            if char != '0' and char != '.':
                start_index = index
                break
        if l[start_index+1] == '.':
            end_index = start_index + 3
        else:
            end_index = start_index + 2
        
        retr = s[start_index:end_index]
        raise_power = (dot_index - start_index - 1) if dot_index - start_index >= 0 else (dot_index - start_index )
        placement = -raise_power # returns position relative to decimal place

        return retr, placement
    
    def get_main_to_placement(num, placement):
        """returns num only up to placement+1 digits"""
        # exception for negative placement longer than non-decimal 
        # component of number
        if placement < -len(str(abs(int(num)))):
            return 0
        
        # body
        f = float(num)
        s = str(f) + '0000000000000'
        l = list(s)
        sdec_beforeoffset = np.asarray(range(len(s)))
        sdec = sdec_beforeoffset - np.argwhere(np.char.strip(l) == '.')[0][0] # position relative to decimal place
        trunc_index = placement + 1 # get the index from placement at where to cut
        trunc_index = -1 if trunc_index == 0 else trunc_index
        end_index = np.argwhere(sdec == trunc_index)[0][0] + 1 # index in num list instead of sdec
        return s[:end_index]

    result = []
    for index, (main_num, error_num) in enumerate(zip(popt, perr)):
        error_str, placement = get_first_two_non_zero_digits(error_num)
        main_str = get_main_to_placement(main_num, placement)
        print(f"B[{index}] = {main_str}({error_str})")
        result.append(f"{main_str}({error_str})")
    return result
    
### dyfit plotting ###
def get_dy(fittype, x, popt, pcov):
    """Returns the dy in fitting for given popt, pcov over x. For use 
    with plt.fill_between.

    Inputs
    ------
    fittype: str
        'linear', 'power2', 'power3'
        denotes the original fitting equation
    x: array of values fitting was done to do filling
    popt, pcov: obtained from curve_fit or ODR_fit
    
    Outputs
    ------
    dy: Array for use in plt.fill_between, such as in
        ax.fill_between(t, f(t, *popt)-dy, f(t, *popt)+dy)
    """
    # convert pcov to to cov to std_dev 
    perr= np.sqrt(np.diag(pcov))

    #dx = 0 for all fitting
    if fittype == 'linear':
        # y = mx+c => dy = dm * x + dc
        dm, dc = perr
        return dm * x + dc
    elif fittype == 'power2':
        # y= a * x**b
        a, b = popt
        da, db = perr
        r = a * (x**b)
        return r* np.sqrt((da/a)**2.0+(db*np.log(x))**2.0)
    elif fittype == 'power3':
        a, b, c = popt
        da, db, dc = perr
        r = a * (x**b)
        return r* np.sqrt((da/a)**2.0+(db*np.log(x))**2.0) + dc
    elif fittype == 'gaussian':
        a, b, c = popt
        da, db, dc = perr
        y = a * np.exp(- (x-b)**2. /c )
        mult = da/a + 2*(x-b)*db/c + (((x-b)/c)**2.)*dc
        return y * mult
    else:
        raise ValueError("fittype specified is unrecognised!")

# save matplotlib figure
def write_fig_to_file(fig, name):
    """Writes matplotlib figure into png.
    Note that an A4 paper in inches is 	8-1/4 x 11-3/4.
    Inputs
    ------
    fig: matplotlib (fig) from fig, axs = plt.subplots(nrows=..., ncols=...)
    name: (string) for filename to be as written as
    """
    DIRECTORY_WRITE = os.path.dirname( __file__ )
    fullpath = os.path.join(DIRECTORY_WRITE, name + '.png')
    fig.savefig(fullpath)
    print(f"figure has been saved to {str(os.path.join(DIRECTORY_WRITE, name + '.png'))}")