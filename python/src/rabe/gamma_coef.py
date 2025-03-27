import numpy as np


def convert_landreman_to_gamma_coefs(
    L31: float,
    L32: float,
    G: float,
    I: float,
    iota: float,
    helicity_n: int,
    dpsi_dr: float,
):

    conversion_fac = (G + helicity_n * I) / (iota - helicity_n) * 1 / 2 / dpsi_dr
    gamma31_redl = conversion_fac * L31
    gamma32_redl = conversion_fac * (L32 + 5 / 2 * L31)
    return gamma31_redl, gamma32_redl


def get_average_gamma_coefs(gamma):
    gamma31 = (gamma[2, 0] + gamma[0, 2]) / 2  # Onsager Symmetry
    gamma31_std = np.sqrt((gamma31 - gamma[2, 0]) ** 2 + (gamma31 - gamma[0, 2]) ** 2)

    gamma32 = (gamma[2, 1] + gamma[1, 2]) / 2  # Onsager Symmetry
    gamma32_std = np.sqrt((gamma32 - gamma[2, 1]) ** 2 + (gamma32 - gamma[1, 2]) ** 2)
    return gamma31, gamma31_std, gamma32, gamma32_std
