# %%
import os
import numpy as np

from simsopt.mhd.vmec import Vmec
from simsopt.mhd.boozer import Boozer
from simsopt.mhd.profiles import ProfilePolynomial
from simsopt.mhd.bootstrap import j_dot_B_Redl, RedlGeomBoozer

import h5py

vmec_file = os.path.join("output/wout_axisymmetric.nc")

neo2_file = os.path.join("output/axisymmetric.out")

with h5py.File(neo2_file, "r") as h5_file:
    boozer_s = h5_file["boozer_s"][()]
    avnabpsi = h5_file["avnabpsi"][()]
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
psi_SI = geom_data.psi_edge

stors = []
gamma31_redl = []
gamma32_redl = []
gamma31_neo2 = []
gamma31_std_neo2 = []
gamma32_neo2 = []
gamma32_std_neo2 = []


# %%
import matplotlib.pyplot as plt

plt.figure()
plt.plot(stors, -np.array(gamma31_redl), "xr", label="Redl, -gamma31")
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
plt.plot(stors, -np.array(gamma32_redl), "xb", label="Redl, -gamma32")
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
