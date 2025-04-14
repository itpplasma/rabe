# %%
import os
import numpy as np

import h5py

from simsopt.mhd.vmec import Vmec
from simsopt.mhd.boozer import Boozer
from simsopt.mhd.bootstrap import RedlGeomBoozer

from rabe.shaing_callen import shaing_callen_bootstrap

magnetic_case = "helical"

if magnetic_case == "axi":
    vmec_file = os.path.join("output/wout_axisymmetric.nc")
    neo2_file = os.path.join("output/axisymmetric_collisionality_scan.out")
    bc_filename = os.path.join("output/axisymmetric.bc")
    helicity_n = 0
elif magnetic_case == "helical":
    vmec_file = os.path.join(
        "output/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc"
    )
    neo2_file = os.path.join(
        "output/helicalsymmetric_collisionality_scan_stor_0p26.out"
    )
    bc_filename = os.path.join("output/quasi_helicalsymmetric.bc")
    helicity_n = -1
    R0 = 1406
else:
    raise ValueError("Unknown magnetic_case")

with h5py.File(neo2_file, "r") as h5_file:
    s_tor_neo2 = h5_file["boozer_s"][()]
    kappa = h5_file["conl_over_mfp"][()]
    avnabs = h5_file["avnabpsi"][()]
    lambda_bB_neo2output = h5_file["alambda_bb"][()]

nu_star = -kappa * np.pi / 2 * R0

eps = 1e-2

vmec = Vmec(vmec_file)
boozer = Boozer(vmec, mpol=16, ntor=16)
geom = RedlGeomBoozer(boozer, [s_tor_neo2[0], s_tor_neo2[0] + eps], helicity_n)

geom_data = geom()
psi_SI_rhs = geom_data.psi_edge
avnabAphi_cgs = avnabs[0] * psi_SI_rhs

iota = geom_data.iota[0]
b_0 = geom_data.Bmax[0]
N = helicity_n * geom_data.nfp

lambda_bB_analytic = shaing_callen_bootstrap(
    bc_filename, s_tor_neo2[0], iota, b_0, avnabAphi_cgs, N=N
)
# %%
import matplotlib.pyplot as plt

plt.figure()
plt.plot(
    nu_star,
    lambda_bB_analytic * np.ones_like(lambda_bB_neo2output),
    "r",
    label="analytic",
)
plt.plot(nu_star, lambda_bB_neo2output, "ob", label="NEO-2")
plt.xlabel("kappa")
plt.ylabel("lambda_bB")
plt.xscale("log")
plt.legend()
plt.show()

plt.figure()
plt.plot(
    nu_star,
    np.array(lambda_bB_analytic) / np.array(lambda_bB_neo2output),
    "or",
    label="analytic/neo2",
)
plt.xlabel("kappa")
plt.ylabel("ratio")
plt.xscale("log")
plt.legend()
plt.show()
