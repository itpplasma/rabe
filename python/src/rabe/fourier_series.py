import numpy as np
from numba import jit, prange
from numba import int64, float64
from numba.experimental import jitclass

fourier_spec = [
    ("m", int64[:]),
    ("n", int64[:]),
    ("fmnc", float64[:]),
    ("fmns", float64[:]),
]


# @jitclass(fourier_spec)
class FourierSeries:
    """
    Class to represent a fourier series with different kernels
    """

    def __init__(self, m, n, fmnc, fmns):
        self.m = m
        self.n = n
        self.fmnc = fmnc
        self.fmns = fmns


# @jit(nopython=True)
def evaluate(fourier: FourierSeries, theta, phi):
    f = np.zeros_like(theta)
    for i in range(0, len(fourier.m)):
        fmnc = fourier.fmnc[i]
        fmns = fourier.fmns[i]
        m = fourier.m[i]
        n = fourier.n[i]
        f += fmnc * np.cos(m * theta + n * phi) + fmns * np.sin(m * theta + n * phi)
    return f
