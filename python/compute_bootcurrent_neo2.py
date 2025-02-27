# %%
import os
import h5py
import numpy as np
from simsopt.mhd.profiles import ProfilePolynomial
from neo2par import get_D_coef_neo2par

s_tor = 0.5

c = 3e10
e = -4.8e-10
EV2ERG = 1.6e-12
GAUSS2TESLA = 1e-4
STATA2AMPERE = 3e-9
STATA_GAUSS2AMPERE_TESLA = STATA2AMPERE * GAUSS2TESLA
CM2_TO_M2 = 1e-4
CM3_TO_M3 = 1e-6

ne_spline = ProfilePolynomial(4.0e20 * CM3_TO_M3 * np.array([1]))
Te_spline = ProfilePolynomial(12.0e3 * EV2ERG * np.array([1, -1]))

neo2par_output = "output"
D_neo2par = get_D_coef_neo2par(neo2par_output, Te_spline(s_tor), m=9.1e-28, Z=-1)

with h5py.File(os.path.join(neo2par_output, "fulltransp.h5")) as output:
    ds_dr = output["avnabpsi"][()]

Te = Te_spline(s_tor)
dTe_ds = Te_spline.dfds(s_tor)
dTe_dr = dTe_ds * ds_dr

ne = ne_spline(s_tor)
dne_ds = ne_spline.dfds(s_tor)
dne_dr = dne_ds * ds_dr

A_1 = 1 / ne * dne_dr - 3 / 2 * 1 / Te * dTe_dr
A_2 = 1 / Te * dTe_dr

j_dot_B = -ne * (D_neo2par[0, 2] * A_1 + D_neo2par[1, 2] * A_2) * e
j_dot_B_SI = j_dot_B * STATA_GAUSS2AMPERE_TESLA / CM2_TO_M2

print("--------------------------------")
print("s = ", s_tor)
print("<j dot B> = ", j_dot_B_SI)
print("<j dot B> = ", j_dot_B_SI * 1e-6, " MegaAmpere Tesla / meter^2")
