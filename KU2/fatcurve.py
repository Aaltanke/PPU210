import numpy as np
def fatcurve(s, c, **kwargs):
    """Extract the cycle amount for a given stress range based of S-N curve in Eurocode 3.

        :param s: Array of stress ranges in MPa.
        :type stress_signal: numpy.ndarray
        :param c: Detail category in MPa.
        :type stress_signal: double

        :return: Array with cycle count for each stress range
        :rtype: numpy.ndarray
        """
    m1 = 3
    m2 = 5 
    d = (2/5) ** (1/m1) * c  # 'd' is the fatigue limit
            
    n = np.zeros(len(s))          
    for i in range(len(n)):
     n1 = c**m1 * 2e6 / s[i]**m1  # assume N <= 5e6
     n2 = d**m2 * 5e6 / s[i]**m2  # assume 5e6 < N <= 1e8
     n3 = float('inf')  # assume N > 1e8
     if n1 <= 5e6:
         n[i] = n1
     elif n2 <= 1e8:
             n[i] = n2
     else:
            n[i] = n3
                    
    return n 