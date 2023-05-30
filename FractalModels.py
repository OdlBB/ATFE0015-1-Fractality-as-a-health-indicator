import numpy as np
from numpy import random
import math


def IvanovModel(N=7, w_SA=0.01, w_SS=0.01, w_PS=0.03, t_SA=0.6, t_SS_lower=0.2, t_SS_higher=1, t_PS_lower=0.9, t_PS_higher=1.5, mean_eta=0, std_eta=0.5*math.sqrt(2), mean_T=1000, std_T=500, Num_beats=20000):
    """
    Simulate a heart rate signal using the Ivanov model

    Arguments:
    - N: the number of sympathetic inputs
    - w_SA: the input feedback strength of the SA node
    - w_SS: the input feedback strength of the sympathetic system(SS)
    - w_PS: the input feedback strength of the parasympathetic system(PS)
    - t_SA: the preferred level of SA
    - t_SS_higher: upper bound of the preferred level of SS
    - t_SS_lower: lower bound of the preferred level of SS
    - t_PS_higher: upper bound of the preferred level of PS
    - t_PS_lower: lower bound of the preferred level of PS
    - mean_eta: the mean of the white noise
    - std_eta: the standard deviation of the white noise
    - mean_T: the mean of the time for a constant preferred level
    - std_T: the standard deviation of the time for a constant preferred level
    - Num_beats: the number of points of the generated signal

    Returns:
    - A numpy array with the generated signal
    """

    # Intialising all variable
    t_n = random.uniform(t_SS_lower, t_SS_higher) #tau(n)
    signal0 = np.zeros(Num_beats) #final signal

    # Time T for a preferred level to remain constant
    T_PS = random.normal(mean_T, std_T)
    T_SS = np.zeros(N)
    for i in range(N):
        T_SS[i] = random.normal(mean_T, std_T)

    cnt_PS = 0 #Used to keep track of whether the prefered level needs to be updated
    cnt_SS = np.zeros(N)

    # Preferred levels
    t_PS = random.uniform(t_PS_lower, t_PS_higher)

    t_SS = np.zeros(N)
    for i in range(N):
        t_SS[i] = random.uniform(t_SS_lower, t_SS_higher)

    
    # Generating the signal
    for n in range(Num_beats):
        diff = 0 # second term of the formula
        
        # SA node
        eta = random.laplace(mean_eta, std_eta/math.sqrt(2))
        if t_n < t_SA:
            diff += w_SA*(1+eta)
        else:
            diff -= w_SA*(1+eta)
            
        # PS
        eta = random.laplace(mean_eta, std_eta/math.sqrt(2))
        
        # Changing the preferred level
        if cnt_PS > T_PS:
            cnt_PS = 0
            T_PS = random.normal(mean_T, std_T)
            t_PS = random.uniform(t_PS_lower, t_PS_higher)
        cnt_PS += 1
        
        if t_n < t_PS:
            diff += w_PS*(1+eta)
        else:
            diff -= w_PS*(1+eta)
            

        # SS
        for j in range(N):
            eta = random.laplace(mean_eta, std_eta/math.sqrt(2))
            
            # Changing the preferred level
            if cnt_SS[j] > T_SS[j]:
                t_SS[j] = random.uniform(t_SS_lower, t_SS_higher)
                cnt_SS[j] = 0
                T_SS[j] = random.normal(mean_T, std_T)
                
            cnt_SS[j] += 1
            
            if t_n < t_SS[j]:
                diff += w_SS*(1+eta)
            else:
                diff -= w_SS*(1+eta)

                
        signal0[n] = t_n
        t_n += diff # Value of tau(n+1)

    return signal0

