# %%
import os
import numpy as np

import h5py

from simsopt.mhd.vmec import Vmec
from simsopt.mhd.boozer import Boozer
from simsopt.mhd.bootstrap import RedlGeomBoozer

from rabe.shaing_callen import shaing_callen_bootstrap


vmec_file = os.path.join("output/wout_axisymmetric.nc")
neo2_file = os.path.join("output/axisymmetric_collisionality_scan.out")
bc_filename = os.path.join("output/axisymmetric.bc")


s_tor_vmec = [0.01, 0.02]
helicity_n = 0
eps = 1e-8

vmec = Vmec(vmec_file)
boozer = Boozer(vmec, mpol=16, ntor=16)
geom = RedlGeomBoozer(boozer, s_tor_vmec, helicity_n)

with h5py.File(neo2_file, "r") as h5_file:
    s_tor_neo2 = h5_file["boozer_s"][()]
    kappa = h5_file["conl_over_mfp"][()]
    avnabs = h5_file["avnabpsi"][()]
    lambda_bB_neo2output = h5_file["alambda_bb"][()]

geom_data = geom()
psi_SI = geom_data.psi_edge
avnabpsi_cgs = avnabs[0] * psi_SI

iota = geom_data.iota[0]
b_0 = geom_data.Bmax[0]

lambda_bB_analytic = shaing_callen_bootstrap(
    bc_filename, s_tor_neo2[0], iota, b_0, avnabpsi_cgs
)

# %%
import matplotlib.pyplot as plt

plt.figure()
plt.plot(
    -kappa,
    -lambda_bB_analytic * np.ones_like(lambda_bB_neo2output),
    "xr",
    label="-analytic",
)
plt.plot(-kappa, lambda_bB_neo2output, "ob", label="NEO-2")
plt.xlabel("kappa")
plt.ylabel("lambda_bB")
plt.xscale("log")
plt.legend()
plt.show()

plt.figure()
plt.plot(
    -kappa,
    np.array(lambda_bB_analytic) / np.array(lambda_bB_neo2output),
    "or",
    label="analytic/neo2",
)
plt.xlabel("kappa")
plt.ylabel("ratio")
plt.xscale("log")
plt.legend()
plt.show()
