"""Tests of the batched VMEC++ adjoint and the least-squares Jacobian (manual, like test_rabe_boundary).

RABE_OPT_INPUT=<simsopt>/tests/test_files/input.LandremanPaul2021_QH_reactorScale_lowres pytest opt/test_lsq_opt.py
"""

import os

import jax
import jax.numpy as jnp
import numpy as np
import pytest
import vmecpp

pytestmark = pytest.mark.skipif("RABE_OPT_INPUT" not in os.environ, reason="set RABE_OPT_INPUT")


def _input(ns=13):
    inp = vmecpp.VmecInput.from_file(os.environ["RABE_OPT_INPUT"])
    return inp.model_copy(update=dict(ns_array=np.array([ns]), ftol_array=np.array([1e-20]),
                                      niter_array=np.array([40000]), always_fix_m1_gauge=True))


def test_batched_adjoint_matches_vmecpp_vjp():
    """vmap over the pullback with one shared factorization equals vmecpp's own VJP per row."""
    import vmecpp_jacobian as vj
    from rabe_boundary import harmonics_of
    inp = _input()
    b0 = jnp.asarray(np.stack([np.asarray(inp.rbc, float), np.asarray(inp.zbs, float)]))

    def plain(b):
        return harmonics_of(vmecpp.run(inp.model_copy(update={"rbc": b[0], "zbs": b[1]}),
                                       max_threads=2, verbose=False).wout)

    def batched(b):
        return harmonics_of(vj.traced_run(inp, b, max_threads=2).wout)

    h, pb = jax.vjp(batched, b0)
    rng = np.random.default_rng(0)
    cots = [{k: jnp.asarray(rng.standard_normal(np.shape(v))) for k, v in h.items()} for _ in range(3)]
    n_fact = vj.STATS.factorizations
    g = vj.batched_pullback(pb, {k: jnp.stack([c[k] for c in cots]) for k in h})
    assert vj.STATS.factorizations == n_fact + 1
    for i, c in enumerate(cots):
        _, pb_ref = jax.vjp(plain, b0)
        ref = np.asarray(pb_ref(c)[0])
        assert np.max(np.abs(g[i] - ref)) <= 1e-9 * np.max(np.abs(ref))


def test_lsq_jacobian_frozen_fd():
    """J d against a central difference with re-solved VMEC++ and the frozen rabe structure."""
    import tempfile

    import lsq_opt as lo
    from alpha_optimize import free_mask, settings_of
    inp = _input()
    st = settings_of(dict(M_pol=-1.0, N_tor=4.0, lines=40, n_gauss=32, n_eta=16, n_gauss_well=16), 4)
    b0 = np.stack([np.asarray(inp.rbc, float), np.asarray(inp.zbs, float)])
    mask = free_mask(b0.shape, int(inp.ntor), 1, 1)
    p = lo.Problem(inp, [0.5], st, 9, 8, "E2", mask, b0, 1.7e-3, lo.Weights(mirror_band=0.01),
                   tempfile.mkdtemp(), max_threads=2)
    rng = np.random.default_rng(1)
    x = 0.3 * rng.standard_normal(p.n)
    ev = p.evaluate(x, "J")
    eg = p.evaluate(x, "g")
    assert np.max(np.abs(eg.g - ev.g)) <= 1e-9 * np.max(np.abs(ev.g))
    d = rng.standard_normal(p.n)
    d /= np.linalg.norm(d)

    def frozen(xx):
        return p.evaluate(xx, structure=ev.structure).r

    h = 0.01
    fd = (frozen(x + h * d) - frozen(x - h * d)) / (2 * h)
    Jd = ev.J @ d
    assert np.max(np.abs(fd - Jd)) <= 1e-3 * np.max(np.abs(Jd))
