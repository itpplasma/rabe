import numpy as np
import os
from simsopt.mhd.vmec import Vmec
from simsopt.mhd.boozer import Boozer, Quasisymmetry

from rabe.boozer_modes import read_modes_bc, split_off_symmetric_modes

output_dir = "output"


def test_split_off_symmetric_modes():
    s_tor_idx = 24
    helicity_n = 1

    bc_file = os.path.join(output_dir, "booz_xform_field.bc")
    get_mode_idx = lambda m, n, n_max: (
        (2 * n_max + 1) * (m - 1) + (n_max - n) + (n_max + 1)
    )
    _, _, _, bmnc = read_modes_bc(bc_file, get_mode_idx)
    _, nonsym_bmnc = split_off_symmetric_modes(bmnc, helicity_n=helicity_n)

    s_tors = nonsym_bmnc.rho_tor**2
    s_tor = s_tors[s_tor_idx]

    vmec_file = os.path.join(
        output_dir, "wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc"
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

    bc_file = os.path.join(output_dir, "booz_xform_field.bc")
    get_mode_idx = lambda m, n, n_max: (
        (2 * n_max + 1) * (m - 1) + (n_max - n) + (n_max + 1)
    )
    _, _, _, bmnc = read_modes_bc(bc_file, get_mode_idx)
    sym_bmnc, _ = split_off_symmetric_modes(bmnc, helicity_n=helicity_n)

    quasi_symmetric_bc = os.path.join(output_dir, "quasi_symmetric.bc")
    _, _, _, quasisym_bmnc = read_modes_bc(quasi_symmetric_bc, get_mode_idx)

    np.testing.assert_allclose(quasisym_bmnc.coefs, sym_bmnc.coefs)


# %%
