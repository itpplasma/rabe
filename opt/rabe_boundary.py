"""Objective and exact gradient of rabe's offset with respect to the VMEC++ boundary.

Chain, per evaluation:

1. VMEC++ fixed-boundary solve of ``(rbc, zbs)`` with the m=1 gauge pinned
   (``vmecpp.run`` with traced boundary, Enzyme build, implicit boundary adjoint,
   reverse mode only).
2. The primary wout arrays ``rmnc``, ``zmns``, ``lmns`` (half mesh), ``iotaf``,
   ``phi[-1]`` and ``Rmajor_p`` of that equilibrium.  The derived ``bmnc``,
   ``bsub*`` and ``gmnc`` never enter.
3. The libneo Boozer transform with the ``-I*lambda`` term (libneo PR 450) in
   JAX, :mod:`pyrabe.boozer_libneo`, which also gives ``|grad s|``, ``psi_edge``
   and ``R``.
4. rabe's derivative mode in JAX, :mod:`pyrabe.spectrum_derivative`, for
   ``Lambda_A``, ``Lambda_B``, ``Lambda_S``, ``nu_star_crit`` and ``eta_b``.

Steps 3 and 4 are :class:`pyrabe.wout_derivative.WoutChain`.  The gradient is
one VMEC++ solve, one VJP of the pyrabe chain and one implicit adjoint solve of
VMEC++.  rabe's discrete choices (``p/q``, field-line labels, bounding maxima)
are frozen at the boundary of the call and recomputed at every new boundary,
so the gradient is that of the smooth branch through the current point.

The default objective is the offset of the 1/nu-regime transport coefficient,
``sum_s w_s (Lambda_A(s) / sqrt(nu_star) + Lambda_B(s) / nu_star)``, as in
``pyrabe.simsopt.RabeOptimizable.Lambda_offset``.  Any other smooth scalar of
the outputs can be supplied as ``reduce`` (a JAX function of the ``(S, 5)``
value array).
"""

from __future__ import annotations

import os
import tempfile
import time
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field

import jax
import jax.numpy as jnp
import numpy as np
import vmecpp
from pyrabe import spectrum_derivative as sd
from pyrabe.wout_derivative import WoutChain

jax.config.update("jax_enable_x64", True)

OUTPUT_NAMES = sd.OUTPUT_NAMES


def offset_reduce(nu_star: float, weights: Sequence[float] | None = None,
                  squared: bool = False) -> Callable:
    """``sum_s w_s L_s`` or ``sum_s w_s L_s^2`` with ``L = Lambda_A / sqrt(nu_star) + Lambda_B / nu_star``."""
    def reduce(values):
        w = jnp.ones(values.shape[0]) if weights is None else jnp.asarray(weights, dtype=float)
        offset = values[:, 0] / jnp.sqrt(nu_star) + values[:, 1] / nu_star
        return jnp.sum(w * (offset**2 if squared else offset))
    return reduce


def harmonics_of(wout) -> dict:
    """Primary arrays of a VMEC++ wout in the layout of :func:`pyrabe.wout_derivative.load_harmonics`."""
    return {
        "rmnc": jnp.asarray(wout.rmnc).T,
        "zmns": jnp.asarray(wout.zmns).T,
        "lmns": jnp.asarray(wout.lmns).T,
        "iotaf": jnp.asarray(wout.iotaf),
        "phiedge": jnp.asarray(wout.phi)[-1],
        "rmajor": jnp.asarray(wout.Rmajor_p),
    }


@dataclass
class Evaluation:
    """Everything one boundary needs: values, the frozen chain and the VMEC++ pullback."""

    boundary: np.ndarray
    harm: dict
    chain: WoutChain
    values: np.ndarray
    pullback: Callable | None
    timing: dict = field(default_factory=dict)


class BoundaryObjective:
    """``objective(boundary)`` and ``grad(boundary)`` for ``boundary = [rbc, zbs]``.

    Args:
        indata: Fixed-boundary, stellarator-symmetric :class:`vmecpp.VmecInput`
            (prescribed iota or current).  Its ``rbc``/``zbs`` are the starting
            boundary; the resolution (``ns_array``, ``ftol_array``) is used as is.
        s_list: Normalized toroidal fluxes of the surfaces.
        settings: rabe settings (helicity ``M_pol``, ``N_tor``, lines, ``n_gauss``).
        mboz, nboz: Boozer resolution, ``m < mboz`` (libneo convention).
        reduce: Scalar of the ``(S, 5)`` outputs; default :func:`offset_reduce`.
        nu_star: Collisionality of the default offset objective.
        squared: Default objective ``sum_s L_s^2`` instead of ``sum_s L_s``.

    ``boundary`` is an array of shape ``(2, mpol, 2 ntor + 1)`` with rows
    ``rbc`` and ``zbs`` in VMEC++'s layout.  The last evaluation is cached, so
    ``objective`` followed by ``grad`` at the same boundary costs one solve.
    """

    def __init__(self, indata, s_list, settings: sd.DerivativeSettings, mboz: int, nboz: int,
                 reduce: Callable | None = None, nu_star: float = 1e-3, squared: bool = False,
                 workdir=None):
        if indata.lasym:
            raise ValueError("stellarator-symmetric equilibria only")
        self.indata = indata.model_copy(update={"always_fix_m1_gauge": True})
        self.s = np.atleast_1d(np.asarray(s_list, dtype=float))
        self.settings = settings
        self.mboz, self.nboz = int(mboz), int(nboz)
        self.reduce = offset_reduce(nu_star, squared=squared) if reduce is None else reduce
        self._dreduce = jax.jit(jax.grad(self.reduce))
        self._workdir = workdir or tempfile.mkdtemp(prefix="rabe_boundary_")
        # The libneo plans depend only on ns, the mode table and nfp, which are the
        # same for every boundary, so one reference wout serves all evaluations.
        self.base_wout = os.path.join(self._workdir, "wout_base.nc")
        vmecpp.run(self.indata, verbose=False).wout.save(self.base_wout)
        self._last: Evaluation | None = None
        self.n_solves = 0

    @property
    def x0(self) -> np.ndarray:
        return np.stack([np.asarray(self.indata.rbc, float), np.asarray(self.indata.zbs, float)])

    def _input(self, boundary):
        return self.indata.model_copy(update={"rbc": boundary[0], "zbs": boundary[1]})

    def _harm_of_boundary(self, boundary):
        return harmonics_of(vmecpp.run(self._input(boundary), verbose=False).wout)

    def evaluate(self, boundary, *, with_pullback: bool = True) -> Evaluation:
        boundary = np.asarray(boundary, dtype=float)
        last = self._last
        if (last is not None and np.array_equal(last.boundary, boundary)
                and (last.pullback is not None or not with_pullback)):
            return last
        t0 = time.perf_counter()
        if with_pullback:
            harm, pullback = jax.vjp(self._harm_of_boundary, jnp.asarray(boundary))
        else:
            harm, pullback = self._harm_of_boundary(boundary), None
        self.n_solves += 1
        t1 = time.perf_counter()
        harm = {k: jnp.asarray(np.asarray(v)) for k, v in harm.items()}
        chain = WoutChain(self.base_wout, self.s, self.settings, self.mboz, self.nboz, harm=harm)
        values = chain.values()
        t2 = time.perf_counter()
        self._last = Evaluation(boundary, harm, chain, values, pullback,
                                {"solve": t1 - t0, "structure_values": t2 - t1})
        return self._last

    def values(self, boundary) -> np.ndarray:
        """``(S, 5)`` outputs with the structure found at this boundary.

        The solve keeps the VMEC++ pullback, so a following :meth:`grad` at the
        same boundary needs only the adjoint solve.
        """
        return self.evaluate(boundary, with_pullback=True).values

    def objective(self, boundary) -> float:
        return float(self.reduce(jnp.asarray(self.values(boundary))))

    def grad(self, boundary) -> np.ndarray:
        """Exact gradient of :meth:`objective`, shape ``(2, mpol, 2 ntor + 1)``."""
        ev = self.evaluate(boundary, with_pullback=True)
        return self.vjp(boundary, np.asarray(self._dreduce(jnp.asarray(ev.values))))

    def vjp(self, boundary, cotangent) -> np.ndarray:
        """``sum_{s,o} cotangent[s, o] d values[s, o] / d boundary``."""
        ev = self.evaluate(boundary, with_pullback=True)
        t0 = time.perf_counter()
        g_harm = ev.chain.vjp(np.asarray(cotangent, dtype=float), ev.harm)
        t1 = time.perf_counter()
        (g,) = ev.pullback({k: jnp.asarray(g_harm[k], dtype=ev.harm[k].dtype) for k in ev.harm})
        t2 = time.perf_counter()
        ev.timing.update({"pyrabe_vjp": t1 - t0, "vmecpp_adjoint": t2 - t1})
        return np.asarray(g)

    def frozen_values(self, boundary, reference: Evaluation) -> np.ndarray:
        """Outputs at ``boundary`` with the discrete structure of ``reference``.

        This is the function whose derivative :meth:`grad` returns, so a finite
        difference of it checks the gradient without discrete switches.
        """
        harm = {k: jnp.asarray(np.asarray(v)) for k, v in self._harm_of_boundary(
            np.asarray(boundary, float)).items()}
        self.n_solves += 1
        return reference.chain.values(harm)
