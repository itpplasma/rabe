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
from rabe_boundary import (AlphaWeights, BoundaryObjective, alpha_reduce,  # noqa: E402
                           aspect_penalty, boundary_geometry, iota_penalty)


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


def test_alpha_objective_with_constraints_matches_resolved_fd():
    """Composite F (all four criteria, normalized) plus aspect and iota penalties."""
    st = sd.DerivativeSettings(M_pol=-1.0, N_tor=4.0, max_n_fieldlines=40, n_gauss=32,
                               n_eta=16, n_gauss_well=16)
    inp = load_input(INPUT, 13)
    nfp, ntor = int(inp.nfp), int(inp.ntor)
    probe = BoundaryObjective(inp, [0.5], st, 7, 6)
    weights = AlphaWeights().normalized(probe.values(probe.x0))
    iota0 = np.asarray(probe._last.harm["iotaf"])
    a0 = float(boundary_geometry(probe.x0, nfp, ntor)[1])
    obj = BoundaryObjective(inp, [0.5], st, 7, 6, reduce=alpha_reduce(weights),
                            constraints=[aspect_penalty(a0 - 0.1, 10.0, nfp, ntor),
                                         iota_penalty(iota0 + 0.01, 100.0)])
    x0 = obj.x0
    F0 = obj.objective(x0)
    assert F0 > 4.0  # four normalized terms plus the active penalties
    g = obj.grad(x0)
    d = smooth_direction(x0, ntor, 7)
    exact = float(np.sum(g * d))
    ref = obj._last
    h = 3e-5

    def total(x):
        v = obj.frozen_values(x, ref)
        return float(obj.reduce(v)) + obj.penalty(x)

    fd = (total(x0 + h * d) - total(x0 - h * d)) / (2 * h)
    assert abs(fd - exact) <= 1e-3 * abs(exact)
