# %%
import os
import numpy as np
from scipy.integrate import quad
import h5py

from simsopt.mhd.vmec import Vmec
from simsopt.mhd.boozer import Boozer
from simsopt.mhd.bootstrap import RedlGeomBoozer

from rabe.fourier_series import FourierSeries, evaluate


def get_bmn_splines(bc_filename):
    from libneo import BoozerFile
    from scipy.interpolate import interp1d

    bc_file = BoozerFile(bc_filename)
    stor = np.array(bc_file.s)
    m = bc_file.m[0]
    n = bc_file.n[0]

    bmnc_spl = interp1d(stor, np.array(bc_file.bmnc), axis=0)
    bmns_spl = interp1d(stor, np.array(bc_file.bmns), axis=0)
    return bmnc_spl, bmns_spl, m, n


def get_Bphi_spline(bc_filename):
    from libneo import BoozerFile
    from scipy.interpolate import interp1d

    c_cgs = 3e10
    bc_file = BoozerFile(bc_filename)
    stor = np.array(bc_file.s)
    Jpol_SI = np.array(bc_file.Jpol_divided_by_nper) * np.array(bc_file.nper)
    Jpol_cgs = Jpol_SI * 3e9
    Bphi_cgs = -2 * Jpol_cgs / c_cgs

    Bphi_spl = interp1d(stor, Bphi_cgs, axis=0)
    return Bphi_spl


vmec_file = os.path.join("output/wout_axisymmetric.nc")
neo2_file = os.path.join("output/axisymmetric.out")
bc_filename = os.path.join("output/axisymmetric.bc")
bmnc_spl, bmns_spl, m, n = get_bmn_splines(bc_filename)
Bphi_spl = get_Bphi_spline(bc_filename)

surfaces = np.linspace(0.013, 0.05, 2)
helicity_n = 0
eps = 1e-8

vmec = Vmec(vmec_file)
boozer = Boozer(vmec, mpol=16, ntor=16)
geom = RedlGeomBoozer(boozer, surfaces, helicity_n)

with h5py.File(neo2_file, "r") as h5_file:
    boozer_s = h5_file["boozer_s"][()]
    avnabs = h5_file["avnabpsi"][()]
    lambda_bB_neo2output = h5_file["alambda_bb"][()]

geom_data = geom()
psi_SI = geom_data.psi_edge

stors = []
lambda_bB_analytic = []
lambda_bB_neo2 = []
for stor in boozer_s:
    if any(np.abs(surfaces - stor) < eps):
        idx_neo2 = np.where(np.abs(boozer_s - stor) < eps)[0][0]

        avnabpsi_cgs = avnabs[idx_neo2] * psi_SI * 1e8

        idx_redl = np.where(np.abs(surfaces - stor) < eps)[0][0]

        iota = geom_data.iota[idx_redl]
        B_0 = geom_data.Bmax[idx_redl]
        bmnc = bmnc_spl(stor).tolist()
        bmns = bmns_spl(stor).tolist()
        fourier_b = FourierSeries(m, n, bmnc, bmns)

        def integrand(x):
            def f(theta):
                B = evaluate(fourier_b, theta, phi=0)
                return np.sqrt(1 - x * B / B_0) / (B / B_0) ** 2

            theta_average, _ = quad(f, 0, 2 * np.pi)
            theta_average = theta_average / (2 * np.pi)
            return x / theta_average

        Bphi_cgs = Bphi_spl(stor)
        integral, _ = quad(integrand, 0, 1)
        lambda_bB = -(3 / 4 * integral - 1) * Bphi_cgs / (avnabpsi_cgs * iota)

        stors.append(stor)
        lambda_bB_analytic.append(lambda_bB)
        lambda_bB_neo2.append(lambda_bB_neo2output[idx_neo2])
# %%
import matplotlib.pyplot as plt

plt.figure()
plt.plot(stors, lambda_bB_analytic, "xr", label="analytic")
plt.plot(stors, lambda_bB_neo2, "ob", label="NEO-2")
plt.xlabel("stor")
plt.ylabel("lambda_bB")
plt.legend()
plt.show()

plt.figure()
plt.plot(
    stors,
    np.array(lambda_bB_analytic) / np.array(lambda_bB_neo2),
    "or",
    label="analytic/neo2",
)
plt.xlabel("stor")
plt.ylabel("ratio")
plt.legend()
plt.show()
