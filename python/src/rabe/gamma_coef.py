import numpy as np


def convert_landreman_to_gamma_coefs(
    G: float, I: float, iota: float, L31: float, L32: float, helicity_n: int
):
    idx_neo2 = np.where(boozer_s == stor)[0][0]

    av = (gamma_out[idx_neo2, 2, 0] + gamma_out[idx_neo2, 0, 2]) / 2
    std = np.sqrt(
        (av - gamma_out[idx_neo2, 2, 0]) ** 2 + (av - gamma_out[idx_neo2, 0, 2]) ** 2
    )
    gamma31_neo2.append(av)
    gamma31_std_neo2.append(std)

    av = (gamma_out[idx_neo2, 2, 1] + gamma_out[idx_neo2, 1, 2]) / 2
    std = np.sqrt(
        (av - gamma_out[idx_neo2, 2, 1]) ** 2 + (av - gamma_out[idx_neo2, 1, 2]) ** 2
    )
    gamma32_neo2.append(av)
    gamma32_std_neo2.append(std)

    dr_ds_cgs = 1 / avnabpsi[idx_neo2]
    dr_ds_SI = dr_ds_cgs * 1e-2
    dr_dpsi_SI = dr_ds_SI / psi_SI

    conversion_fac = (G + helicity_n * I) / (iota - helicity_n) * 1 / 2 * dr_dpsi_SI

    gamma31_redl = conversion_fac * L31
    gamma32_redl = conversion_fac * (L32 + 5 / 2 * L31)
