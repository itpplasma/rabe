# %%
import os
import numpy as np

from simsopt.mhd.vmec import Vmec
from simsopt.mhd.boozer import Boozer
from simsopt.mhd.profiles import ProfilePolynomial
from simsopt.mhd.bootstrap import j_dot_B_Redl, RedlGeomBoozer

filename = os.path.join(
    "python/output/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc"
)

ne = ProfilePolynomial(4.0e20 * np.array([1]))
Te = ProfilePolynomial(12.0e3 * np.array([1, -1]))
Ti = ProfilePolynomial(6.0e3 * np.array([1]))
surfaces = np.linspace(0.3, 0.7, 3)
Zeff = 1
helicity_n = -1

vmec = Vmec(filename)
boozer = Boozer(vmec, mpol=16, ntor=16)
geom = RedlGeomBoozer(boozer, surfaces, helicity_n)
jdotB, details = j_dot_B_Redl(ne, Te, Ti, Zeff, geom=geom, plot=False)

geom_data = geom()
G = geom_data.G
I = geom_data.I
iota = geom_data.iota

# for stor = 0.5
dr_ds_cgs = 1 / 0.00774091
dr_ds_SI = dr_ds_cgs * 1e-2
psi_SI = 4.186388e01 / 2 * np.pi
dr_dpsi_SI = dr_ds_SI / psi_SI

conversion_fac = -(G + helicity_n * I) / (iota - helicity_n) * 1 / 2 * dr_dpsi_SI
gamma_31 = conversion_fac * details.L31
gamma_32 = conversion_fac * (details.L32 + 5 / 2 * details.L31)

import matplotlib.pyplot as plt

plt.scatter(surfaces, jdotB, label="Redl, f_t from Boozer")
plt.xlabel("s")
plt.title("J dot B")
plt.legend(loc=0)
plt.show()

for idx in range(len(surfaces)):
    print("--------------------------------")
    print("s = ", surfaces[idx])
    print("L31 = ", details.L31[idx])
    print("L32 = ", details.L32[idx])
    print("L34 = ", details.L34[idx])
    print("alpha = ", details.alpha[idx])
    print("gamma_31 = ", gamma_31[idx])
    print("gamma_32 = ", gamma_32[idx])
    print("<j dot B> = ", jdotB[idx])
