import os
import copy
import numpy as np
import matplotlib.pyplot as plt
from simsopt.mhd.vmec import Vmec
from simsopt.mhd.boozer import Boozer, Quasisymmetry

from rabe.boozer_modes import read_modes_bc, split_off_symmetric_modes
from rabe.boozer_modes import get_xyz_surface, Modes


OUTPUT_DIR = "output"


def test_split_off_symmetric_modes():
    s_tor_idx = 24
    helicity_n = 1

    bc_file = os.path.join(OUTPUT_DIR, "booz_xform_field.bc")
    get_mode_idx = lambda m, n, n_max: (
        (2 * n_max + 1) * (m - 1) + (n_max - n) + (n_max + 1)
    )
    _, _, _, bmnc = read_modes_bc(bc_file, get_mode_idx)
    _, nonsym_bmnc = split_off_symmetric_modes(bmnc, helicity_n=helicity_n)

    s_tors = nonsym_bmnc.rho_tor**2
    s_tor = s_tors[s_tor_idx]

    vmec_file = os.path.join(
        OUTPUT_DIR, "wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc"
    )
    vmec = Vmec(vmec_file)
    boozer = Boozer(vmec)
    boozer.register(s_tor)
    boozer.run()

    error = Quasisymmetry(boozer, s_tor, helicity_m=1, helicity_n=helicity_n).J()

    np.testing.assert_allclose(
        sum(np.abs(error)),
        np.sum(np.abs(nonsym_bmnc.coefs[s_tor_idx, :] / bmnc.coefs[s_tor_idx, 0])),
    )


def test_quasi_symmetric_bc():
    helicity_n = 1

    bc_file = os.path.join(OUTPUT_DIR, "booz_xform_field.bc")
    get_mode_idx = lambda m, n, n_max: (
        (2 * n_max + 1) * (m - 1) + (n_max - n) + (n_max + 1)
    )
    _, _, _, bmnc = read_modes_bc(bc_file, get_mode_idx)
    sym_bmnc, _ = split_off_symmetric_modes(bmnc, helicity_n=helicity_n)

    quasi_symmetric_bc = os.path.join(OUTPUT_DIR, "quasi_symmetric.bc")
    _, _, _, quasisym_bmnc = read_modes_bc(quasi_symmetric_bc, get_mode_idx)

    np.testing.assert_allclose(quasisym_bmnc.coefs, sym_bmnc.coefs)


def test_get_xyz_surface():
    rmnc = Modes()
    rmnc.coefs = np.array([[1.0, 0.5]])
    rmnc.m = np.array([[0, 1]])
    rmnc.n = np.array([[0, 0]])
    rmnc.get_mode_idx = lambda m, n: m + n

    zmns = copy.deepcopy(rmnc)
    zmns.coefs = np.array([[0.0, 0.5]])
    vmns = copy.deepcopy(rmnc)
    vmns.coefs = np.array([[0.0, 0.0]])
    bmnc = copy.deepcopy(rmnc)
    bmnc.coefs = np.array([[1.0, -1.0]])

    x, y, z, B = get_xyz_surface(rmnc, zmns, vmns, bmnc, 1, n_phi=10, n_theta=10)
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")
    norm = plt.Normalize(B.min(), B.max())
    colors = plt.cm.viridis(norm(B))
    ax.plot_surface(x, y, z, facecolors=colors, edgecolor="none")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.axis("equal")
    mappable = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
    mappable.set_array(B)
    fig.colorbar(mappable, ax=ax, shrink=0.5, aspect=10)
    plt.show()
