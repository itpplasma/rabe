"""Boundary gradient against a re-solved central difference (needs vmecpp with Enzyme and pyrabe).

Run with ``RABE_OPT_INPUT=<input.LandremanPaul2021_QH_reactorScale_lowres> pytest opt``.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

vmecpp = pytest.importorskip("vmecpp")
pytest.importorskip("pyrabe.wout_derivative")
if not vmecpp.has_exact_force_jacobian():
    pytest.skip("vmecpp without the Enzyme force Jacobian", allow_module_level=True)
INPUT = os.environ.get("RABE_OPT_INPUT")
if not INPUT:
    pytest.skip("set RABE_OPT_INPUT to the LP QH reactor-scale lowres input", allow_module_level=True)

from pyrabe import spectrum_derivative as sd  # noqa: E402

from check_gradient import load_input, smooth_direction  # noqa: E402
from rabe_boundary import BoundaryObjective  # noqa: E402


def test_boundary_gradient_matches_resolved_fd():
    st = sd.DerivativeSettings(M_pol=-1.0, N_tor=4.0, max_n_fieldlines=40, n_gauss=32)
    inp = load_input(INPUT, 13)
    obj = BoundaryObjective(inp, [0.5], st, 7, 6, nu_star=1e-3)
    x0 = obj.x0
    obj.objective(x0)
    g = obj.grad(x0)
    assert obj.n_solves == 1
    d = smooth_direction(x0, int(inp.ntor), 7)
    exact = float(np.sum(g * d))
    ref = obj._last
    h = 3e-5
    fd = float((obj.reduce(obj.frozen_values(x0 + h * d, ref))
                - obj.reduce(obj.frozen_values(x0 - h * d, ref))) / (2 * h))
    assert abs(fd - exact) <= 1e-3 * abs(exact)
