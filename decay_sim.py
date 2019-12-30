
import numpy as np
import scipy as sp
from scipy import constants as c
from numba import njit
from numba import prange

@njit
def do_sum(period=9.66, n=int(10**13), r=10**-12):
    total = 0
    for i in range(n):
        total += period * (1 - r)**(3*n/2)
    return total

@njit(fastmath=True)
def do_harder_sum(n=int(6.1525*10**13), dE=10**7, m_hippo=10**13, r0=385000000):
    total = 0
    G = c.G
    m_earth = 5.972 * 10**24
    r = r0
    n = 1
    while r >= 6371000:
        n += 1
        r = (G * m_earth * m_hippo * r0) / (dE * n * r0 + G * m_earth * m_hippo)
        # divide by two because two encounters per orbit, then conver to years
        total += np.sqrt((4*np.pi**2 / (G * m_earth)) * r**3) / (2 * 60*60*24*365)
    print(r)
    print(n)
    return total


   
