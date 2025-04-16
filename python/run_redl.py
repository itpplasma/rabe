# %%
import os
import numpy as np

from simsopt.mhd.vmec import Vmec
from simsopt.mhd.boozer import Boozer
from simsopt.mhd.profiles import ProfilePolynomial
from simsopt.mhd.bootstrap import j_dot_B_Redl, RedlGeomBoozer

filename = os.path.join(
    "output/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc"
)

ne = ProfilePolynomial(4.13e20 * np.array([1]))
Te = ProfilePolynomial(12.0e3 * np.array([1, -1]))
Ti = Te
surfaces = np.linspace(0.1, 0.9, 30)
Zeff = 1
helicity_n = -1

vmec = Vmec(filename)
boozer = Boozer(vmec, mpol=16, ntor=16)
geom = RedlGeomBoozer(boozer, surfaces, helicity_n)
jdotB, details = j_dot_B_Redl(ne, Te, Ti, Zeff, geom=geom, plot=False)

import matplotlib.pyplot as plt

plt.scatter(surfaces, jdotB, label="Redl, f_t from Boozer")
plt.xlabel("s")
plt.title("J dot B")
plt.legend(loc=0)
plt.show()

if False:
    for idx in range(len(surfaces)):
        print("--------------------------------")
        print("s = ", surfaces[idx])
        print("L31 = ", details.L31[idx])
        print("L32 = ", details.L32[idx])
        print("L34 = ", details.L34[idx])
        print("alpha = ", details.alpha[idx])
        print("<j dot B> = ", jdotB[idx])
