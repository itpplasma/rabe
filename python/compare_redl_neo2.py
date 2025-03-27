# %%
import os
import numpy as np

from simsopt.mhd.vmec import Vmec
from simsopt.mhd.boozer import Boozer
from simsopt.mhd.profiles import ProfilePolynomial
from simsopt.mhd.bootstrap import j_dot_B_Redl, RedlGeomBoozer

import h5py

from rabe.gamma_coef import get_average_gamma_coefs, convert_landreman_to_gamma_coefs

vmec_file = os.path.join("output/wout_axisymmetric.nc")

neo2_file = os.path.join("output/axisymmetric.out")

with h5py.File(neo2_file, "r") as h5_file:
    boozer_s = h5_file["boozer_s"][()]
    ds_dr = h5_file["avnabpsi"][()]
    gamma_out = h5_file["gamma_out"][()]


ne = ProfilePolynomial(4.0e20 * np.array([1]))
Te = ProfilePolynomial(12.0e3 * np.array([1, -1]))
Ti = ProfilePolynomial(6.0e3 * np.array([1]))

Zeff = 1
helicity_n = -0

vmec = Vmec(vmec_file)
boozer = Boozer(vmec, mpol=16, ntor=16)
geom = RedlGeomBoozer(boozer, boozer_s, helicity_n)
jdotB, details = j_dot_B_Redl(ne, Te, Ti, Zeff, geom=geom, plot=False)

geom_data = geom()
L31 = details.L31
L32 = details.L32
G = geom_data.G
I = geom_data.I
iota = geom_data.iota
dpsi_dr = geom_data.psi_edge * ds_dr

gamma31_redl = np.zeros_like(boozer_s)
gamma32_redl = np.zeros_like(boozer_s)
gamma31_neo2 = np.zeros_like(boozer_s)
gamma31_std_neo2 = np.zeros_like(boozer_s)
gamma32_neo2 = np.zeros_like(boozer_s)
gamma32_std_neo2 = np.zeros_like(boozer_s)

for idx in range(boozer_s):
    (
        gamma31_neo2[idx],
        gamma31_std_neo2[idx],
        gamma32_neo2[idx],
        gamma32_std_neo2[idx],
    ) = get_average_gamma_coefs(gamma_out[idx])
    gamma31_redl[idx], gamma32_redl[idx] = convert_landreman_to_gamma_coefs(
        L31[idx], L32[idx], G[idx], I[idx], iota[idx], helicity_n, dpsi_dr[idx]
    )


# %%
import matplotlib.pyplot as plt

plt.figure()
plt.plot(boozer_s, -np.array(gamma31_redl), "xr", label="Redl, -gamma31")
plt.errorbar(
    boozer_s,
    gamma31_neo2,
    gamma31_std_neo2,
    fmt="o",
    markerfacecolor="red",
    markeredgecolor="none",
    ecolor="red",
    label="Neo2, gamma31",
)
plt.plot(boozer_s, -np.array(gamma32_redl), "xb", label="Redl, -gamma32")
plt.errorbar(
    boozer_s,
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
plt.plot(
    boozer_s, np.array(gamma31_redl) / np.array(gamma31_neo2), "or", label="gamma31"
)
plt.plot(
    boozer_s, np.array(gamma32_redl) / np.array(gamma32_neo2), "ob", label="gamma32"
)
plt.xlabel("stor [1]")
plt.ylabel("ratio [1]")
plt.title("comparison dimensional bootstrap coefs")
plt.legend(loc=0)
plt.show()
