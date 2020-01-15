#!/usr/bin/python

import numpy as np
from matplotlib import pyplot as plt


def main(r=.1, R=1, w=1/3300, n=100):
    # out total time will be n periods
    T = n/w 
    # lets just do 100 points per period
    t = np.arange(0, T, .01/w)
    theta = 2 * np.pi * w * t
    # calculate distances
    d = dist(theta, r, R)

    # now take fourier transform
    fft_d = np.fft.rfft(d)
    freqs = np.fft.rfftfreq(len(t), 0.01/w)
    plt.plot(freqs, fft_d.real)
    plt.plot(freqs, fft_d.imag)
    return fft_d, freqs


def dist(theta, r=.1, R=1):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    # let the observer be at x=R, y=0)
    return np.sqrt((R - x) ** 2 + y ** 2)
