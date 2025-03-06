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


@jitclass(fourier_spec)
class FourierSeries:
    """
    Class to represent a fourier series
    """

    def __init__(self, m: list, n: list, fmnc: list, fmns: list):
        self.m = np.array(m, dtype=np.int64)
        self.n = np.array(n, dtype=np.int64)
        self.fmnc = np.array(fmnc, dtype=np.float64)
        self.fmns = np.array(fmns, dtype=np.float64)


@jit(nopython=True)
def evaluate(fourier: FourierSeries, theta, phi):
    if not (np.shape(theta) == np.shape(phi)):
        ValueError("theta and phi need to be the same shape")
    f = np.zeros_like(theta)
    for i, _ in enumerate(fourier.m):
        fmnc = fourier.fmnc[i]
        fmns = fourier.fmns[i]
        m = fourier.m[i]
        n = fourier.n[i]
        f += fmnc * np.cos(m * theta - n * phi) + fmns * np.sin(m * theta - n * phi)
    return f
