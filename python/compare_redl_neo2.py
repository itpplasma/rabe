# %%
import os
import numpy as np

from simsopt.mhd.vmec import Vmec
from simsopt.mhd.boozer import Boozer
from simsopt.mhd.profiles import ProfilePolynomial
from simsopt.mhd.bootstrap import j_dot_B_Redl, RedlGeomBoozer

import h5py

vmec_file = os.path.join(
    "output/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc"
)

neo2_file = os.path.join("output/quasi_symmetric.out")

psi_SI = 4.186388e01 / 2 * np.pi  # toroidal magnetic flux at separatrix

ne = ProfilePolynomial(4.0e20 * np.array([1]))
Te = ProfilePolynomial(12.0e3 * np.array([1, -1]))
Ti = ProfilePolynomial(6.0e3 * np.array([1]))
surfaces = np.linspace(0.1, 0.9, 17)
Zeff = 1
helicity_n = -1

vmec = Vmec(vmec_file)
boozer = Boozer(vmec, mpol=16, ntor=16)
geom = RedlGeomBoozer(boozer, surfaces, helicity_n)
jdotB, details = j_dot_B_Redl(ne, Te, Ti, Zeff, geom=geom, plot=False)

with h5py.File(neo2_file, "r") as h5_file:
    boozer_s = h5_file["boozer_s"][()]
    avnabpsi = h5_file["avnabpsi"][()]
    gamma_out = h5_file["gamma_out"][()]

geom_data = geom()

stors = []
gamma31_redl = []
gamma32_redl = []
gamma31_neo2 = []
gamma31_std_neo2 = []
gamma32_neo2 = []
gamma32_std_neo2 = []
for stor in boozer_s:
    if stor in surfaces:
        idx_neo2 = np.where(boozer_s == stor)[0][0]

        av = (gamma_out[idx_neo2, 2, 0] + gamma_out[idx_neo2, 0, 2]) / 2
        std = np.sqrt(
            (av - gamma_out[idx_neo2, 2, 0]) ** 2
            + (av - gamma_out[idx_neo2, 0, 2]) ** 2
        )
        gamma31_neo2.append(av)
        gamma31_std_neo2.append(std)

        av = (gamma_out[idx_neo2, 2, 1] + gamma_out[idx_neo2, 1, 2]) / 2
        std = np.sqrt(
            (av - gamma_out[idx_neo2, 2, 1]) ** 2
            + (av - gamma_out[idx_neo2, 1, 2]) ** 2
        )
        gamma32_neo2.append(av)
        gamma32_std_neo2.append(std)

        dr_ds_cgs = 1 / avnabpsi[idx_neo2]
        dr_ds_SI = dr_ds_cgs * 1e-2
        dr_dpsi_SI = dr_ds_SI / psi_SI

        idx_redl = np.where(surfaces == stor)[0][0]

        G = geom_data.G[idx_redl]
        I = geom_data.I[idx_redl]
        iota = geom_data.iota[idx_redl]
        L31 = details.L31[idx_redl]
        L32 = details.L32[idx_redl]

        conversion_fac = (G + helicity_n * I) / (iota - helicity_n) * 1 / 2 * dr_dpsi_SI

        stors.append(stor)
        gamma31_redl.append(conversion_fac * L31)
        gamma32_redl.append(conversion_fac * (L32 + 5 / 2 * L31))


import matplotlib.pyplot as plt

plt.figure()
plt.plot(stors, gamma31_redl, "xr", label="Redl, gamma31")
plt.errorbar(
    stors,
    gamma31_neo2,
    gamma31_std_neo2,
    fmt="o",
    markerfacecolor="red",
    markeredgecolor="none",
    ecolor="red",
    label="Neo2, gamma31",
)
plt.plot(stors, gamma32_redl, "xb", label="Redl, gamma32")
plt.errorbar(
    stors,
    gamma32_neo2,
    gamma32_std_neo2,
    fmt="o",
    markerfacecolor="blue",
    markeredgecolor="none",
    ecolor="blue",
    label="Neo2, gamma32",
)
plt.xlabel("stor [1]")
plt.ylabel("gamma [1]")
plt.title("comparison dimensional bootstrap coefs")
plt.legend(loc=0)
plt.show()

plt.figure()
plt.plot(stors, np.array(gamma31_redl) / np.array(gamma31_neo2), "or", label="gamma31")
plt.plot(stors, np.array(gamma32_redl) / np.array(gamma32_neo2), "ob", label="gamma32")
plt.xlabel("stor [1]")
plt.ylabel("ratio [1]")
plt.title("comparison dimensional bootstrap coefs")
plt.legend(loc=0)
plt.show()
