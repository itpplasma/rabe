# %%
import sys
import os
import h5py
import numpy as np
from simsopt.mhd.profiles import ProfilePolynomial

s_tor = 0.5
c = 3e10
plasma_radius = 170
EV2ERG = 1.6e-12
GAUSS2TESLA = 1e-4
STATA2AMPERE = 3e-9
STATA_GAUSS2AMPERE_TESLA = STATA2AMPERE * GAUSS2TESLA
CM2_TO_M2 = 1e-4
CM3_TO_M3 = 1e-6

if len(sys.argv) >= 2:
    neo2_ouput = os.path.join(sys.argv[1], "efinal.h5")
    neo2_transport = os.path.join(sys.argv[1], "fulltransp.h5")
else:
    neo2_ouput = "output/efinal.h5"
    neo2_transport = "output/fulltransp.h5"

neo2_ouput = "output/efinal.h5_lorentz"
neo2_transport = "output/fulltransp.h5_lorentz"

with h5py.File(neo2_ouput) as output:
    boot_coef = output["alambda_bb"][()]
with h5py.File(neo2_transport) as output:
    dr_ds = 1 / output["avnabpsi"][()]

ne_spline = ProfilePolynomial(4.0e20 * CM3_TO_M3 * np.array([1]))
Te_spline = ProfilePolynomial(12.0e3 * np.array([1, -1]))

dTe_ds = Te_spline.dfds(s_tor)
dTe_dr = dTe_ds * 1 / dr_ds * EV2ERG

ne = ne_spline(s_tor)

j_dot_B = -c * boot_coef * ne * dTe_dr
j_dot_B_SI = j_dot_B * STATA_GAUSS2AMPERE_TESLA / CM2_TO_M2

print("--------------------------------")
print("s = ", s_tor)
print("<j dot B> = ", j_dot_B_SI)
print("<j dot B> = ", j_dot_B_SI * 1e-6, " MegaAmpere Tesla / meter^2")
