import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import signal
from MFDFA import MFDFA

"""
First method for DFA
"""

def step1(s):
    """
    Convert a signal into a mean-centered cumulative sum.

    Arguments:
    - s: the signal to be converted

    Returns:
    - the converted signal
    """
    # Remove the mean
    mean = np.mean(s)
    mean_array = np.ones(len(s))*mean
    s = np.subtract(s, mean_array)

    # Compute the cumumative sums
    cumulative_sum = np.zeros(len(s))

    tmp = 0
    
    for i in range(len(s)):
        tmp += s[i]
        cumulative_sum[i] = tmp
    
    return cumulative_sum


def rms(s, scale):
    """
    Divides a signal into smaller sub-signals of the same size, compute the 
    root-mean squared of each sub-signal detrended and average all the results.

    Arguments:
    - s: the signal on which the operations are made
    - scale: the size of the sub-signals

    Returns:
    - the mean of the RMS of each sub-signals
    """
    
    new_scale = int(scale)
    
    div_signal = []
    start = 0
    end = len(s)
    
    # Divide the signal into sub-signals of size scale
    for i in range(start, end, new_scale):
        if i+new_scale < end:
            div_signal.append(s[i:i+new_scale])
    
    
    rmss = np.zeros(len(div_signal))
    i = 0
    
    # Compute the RMS of each sub-signal
    for x in div_signal:
        detrended = signal.detrend(x)
        rmss[i] = math.sqrt(np.mean(np.square(detrended)))
        i = i + 1
        
    return np.mean(rmss)
        
    
def DFA(s, plot=False, limit_end=500, limit_start=5, num_points=30):
    """
    Compute the alpha coefficient of a signal using the DFA method.

    Arguments:
    - s: the signal on which the alpha coefficient is computed
    - plot: True to plot the curve of the DFA method
    - limit_end: the largest scale to be considered
    - limit_start: the smallest scale to be considered
    - num_points: the number of scale to consider (and thus the number of 
                    points on which the alpha coefficient is computed)

    Returns:
    - the alpha coefficient
    """
    
    if len(s)<limit_start or limit_start <= 0:
        return 0
    
    if len(s) < limit_end:
        limit_end = int(len(s)/8)

    # Define thte scales
    x = np.linspace(np.log(limit_start), np.log(limit_end), num_points)
    scales = np.exp(x)

    s = step1(s) # Transform the signal into a cumulativ sum
    rmss = np.zeros(len(scales)) 
    i = 0

    # Compute the mean RMS of all scales
    for scale in scales:
        rmss[i] = rms(s, scale) 
        i = i + 1
        
    # Put it on log format
    log_rmss = np.log(rmss)
    log_scales = np.log(scales)
    
    # Compute alpha
    coefs = np.polyfit(log_scales, log_rmss, deg=1)

    if plot==True:
        y = x*coefs[0] + coefs[1]
        
        plt.plot(log_scales, log_rmss, 'o')
        plt.plot(x, y)
        plt.plot()
        plt.xlabel('scales(log)')
        plt.ylabel(r'$\alpha$(log)')
        plt.show()
    
    return coefs[0]


"""
First method but with different tools to speed up the computation
"""

def detrend(sig, order):
    """
    Detrend a signal

    Arguments:
    - sig: the signal to be detrended
    - order: the order of the detrending

    Returns:
    - the detrended signal
    """
    coef = np.polyfit(range(len(sig)), sig, order)
    detrended_signal = sig - np.polyval(coef, range(len(sig)))

    return detrended_signal


def rms_fast(s, scale, order):
    """
    Divides a signal into smaller sub-signals of the same size, compute the 
    root-mean squared of each sub-signal detrended and average all the results.

    Arguments:
    - s: the signal on which the operations are made
    - scale: the size of the sub-signals
    - order: the order of the detrending

    Returns:
    - the mean of the RMS of each sub-signals
    """

    new_scale = int(scale)

    # Divide the signal into sub-signals of size scale
    div_signal = np.array([s[i-new_scale:i] for i in range(new_scale, len(s), new_scale)])
    
    rmss = np.zeros(len(div_signal))
    i = 0
    
    # Compute the RMS of each sub-signal
    for x in div_signal:
        detrended = detrend(x, order)
        rmss[i] = math.sqrt(np.mean(np.square(detrended)))
        i = i + 1
        
    return np.mean(rmss)
        
    
def DFA_fast(s, plot=False, limit_end=500, limit_start=5, num_points=30, order=1):
    """
    Compute the alpha coefficient of a signal using the DFA method.

    Arguments:
    - s: the signal on which the alpha coefficient is computed
    - plot: True to plot the curve of the DFA method
    - limit_end: the largest scale to be considered
    - limit_start: the smallest scale to be considered
    - num_points: the number of scale to consider (and thus the number of 
                    points on which the alpha coefficient is computed)

    Returns:
    - the alpha coefficient
    """
    
    if len(s)<limit_start or limit_start <= 0:
        return 0
    
    if len(s) < limit_end:
        limit_end = int(len(s)/8)
        

    # Define the scales
    x = np.linspace(np.log(limit_start), np.log(limit_end), num_points)
    scales = np.exp(x)

    s = np.cumsum(s - np.mean(s)) # Convert the signal into a mean-centered cumulative sum
    rmss = np.zeros(len(scales))
    i = 0
    
    # Compute the mean RMS of each scale
    for scale in scales:
        rmss[i] = rms_fast(s, scale, order)
        i = i + 1
        
    log_rmss = np.log(rmss)
    log_scales = np.log(scales)
    
    # Compute alpha
    coefs = np.polyfit(log_scales, log_rmss, deg=1)

    if plot==True:
        y = x*coefs[0] + coefs[1]
        
        plt.plot(log_scales, log_rmss, 'o')
        plt.plot(x, y)
        plt.plot()
        plt.xlabel('scales(log)')
        plt.ylabel(r'$\alpha$(log)')
        #plt.title('DFA analysis of subject ') #to change
        plt.show()
    
    return coefs[0]


"""
Second method for DFA
"""

def DFA2(s, lag, plot=False, title="", order=1, save=False, show=True):
    """
    Compute the alpha coefficient of a signal using the DFA method.

    Arguments:
    - s: the signal on which the alpha coefficient is computed
    - lag: a list of the scales considered
    - plot: True to plot the curve of the DFA method
    - title: the title of the plot
    - order: the order of the detrend
    - save: True to save the plot to the pdf format
    - show: True to show the plot

    Returns:
    - the alpha coefficient
    """
    q = 2

    lag2, dfa = MFDFA(np.array(s), lag = lag, q = q, order = order)
    H_hat = np.polyfit(np.log(lag2)[4:30],np.log(dfa[4:30]),1)
    alpha2 = H_hat[0][0]

    # Plot
    y = []
    log_rmss = np.log(dfa)
    log_scales = np.log(lag2)

    for i in range(len(log_scales)):
        y.append(log_scales[i] * H_hat[0][0] + H_hat[1][0])

    if plot:
        plt.plot(log_scales, y)
        plt.plot(log_scales, log_rmss, 'o')
        plt.xlabel("scale(log)")
        plt.ylabel(r'$\alpha$(log)')
        plt.title(title)
        if save:
            plt.savefig(title + ".pdf")
        if show:
            plt.show()

    return alpha2


"""
Power spectrum analysis
"""

def PowerSpectrumAnalysis(data1, plot=False):
    """
    Compute the beta coefficient of the power spectrum analysis method

    Arguments:
    - data1: the signal on which the beta coefficient is computed
    - plot: True to plot the curve

    Returns:
    - the beta coefficient
    """
    time_step = 1 / 1000
    ps1 = np.abs(np.fft.fft(data1)) ** 2
    freqs1 = np.fft.fftfreq(len(data1), time_step)

    # remove the DC bin because 1/f cannot fit at 0 Hz
    ps1 = ps1[freqs1 != 0]
    freqs1 = freqs1[freqs1 != 0]

    idx1 = np.argsort(freqs1)

    f1 = freqs1[idx1]
    p1 = ps1[idx1]

    f2 = np.log10(f1[int(len(f1)/2)+1:])
    p3 = np.log10(p1[int(len(p1)/2)+1:])


    # Compute beta
    coefs = np.polyfit(f2, p3, deg=1)

    # Plot
    if plot==True:
        x = np.linspace(0,5,100)
        y = x*coefs[0] + coefs[1]

        plt.plot(f2, p3)
        plt.plot(x, y)
        plt.xlim(0,3)
        plt.xlabel('f')
        plt.ylabel('S(f)')
        plt.show()

    
    return -coefs[0]