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
   ``Lambda_A``, ``Lambda_B``, ``Lambda_S``, ``nu_star_crit``, ``eta_b`` and
   the alpha criteria ``Lambda_alpha`` and ``Gamma_well``.

Steps 3 and 4 are :class:`pyrabe.wout_derivative.WoutChain`.  The gradient is
one VMEC++ solve, one VJP of the pyrabe chain and one implicit adjoint solve of
VMEC++.  rabe's discrete choices (``p/q``, field-line labels, bounding maxima)
are frozen at the boundary of the call and recomputed at every new boundary,
so the gradient is that of the smooth branch through the current point.

The default objective is the offset of the 1/nu-regime transport coefficient,
``sum_s w_s (Lambda_A(s) / sqrt(nu_star) + Lambda_B(s) / nu_star)``, as in
``pyrabe.simsopt.RabeOptimizable.Lambda_offset``.  Any other smooth scalar of
the outputs can be supplied as ``reduce`` (a JAX function of the ``(S, 7)``
value array), for example the alpha objective :func:`alpha_reduce`.
``constraints`` adds penalties that are JAX functions of the boundary and of
the primary wout arrays (:func:`aspect_penalty`, :func:`volume_penalty`,
:func:`iota_penalty`).
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
N_OUT = len(OUTPUT_NAMES)
I_OUT = {name: i for i, name in enumerate(OUTPUT_NAMES)}


def offset_reduce(nu_star: float, weights: Sequence[float] | None = None,
                  squared: bool = False) -> Callable:
    """``sum_s w_s L_s`` or ``sum_s w_s L_s^2`` with ``L = Lambda_A / sqrt(nu_star) + Lambda_B / nu_star``."""
    def reduce(values):
        w = jnp.ones(values.shape[0]) if weights is None else jnp.asarray(weights, dtype=float)
        offset = values[:, 0] / jnp.sqrt(nu_star) + values[:, 1] / nu_star
        return jnp.sum(w * (offset**2 if squared else offset))
    return reduce


@dataclass(frozen=True)
class AlphaWeights:
    """Weights of the alpha objective :func:`alpha_reduce` (Kasilov's criteria).

    ``nu_star_crit`` aligns the maxima (1), ``offset`` weights ``Lambda_A^2 +
    Lambda_B^2`` (2a), ``Lambda_alpha`` its square (2b) and ``Gamma_well`` the
    drift within the wells (3).  ``reference`` divides each term by its value
    at a reference boundary (``AlphaWeights.normalized``), so that the weights
    are relative.
    """

    nu_star_crit: float = 1.0
    offset: float = 1.0
    Lambda_alpha: float = 1.0
    Gamma_well: float = 1.0
    reference: tuple[float, float, float, float] = (1.0, 1.0, 1.0, 1.0)

    def normalized(self, values) -> "AlphaWeights":
        """Same weights with the reference set to the terms at ``values`` ``(S, 7)``."""
        terms = np.asarray(alpha_terms(jnp.asarray(values)))
        ref = tuple(float(t) if t > 0.0 else 1.0 for t in terms.sum(axis=0))
        return AlphaWeights(self.nu_star_crit, self.offset, self.Lambda_alpha,
                            self.Gamma_well, ref)


def alpha_terms(values):
    """``(S, 4)`` terms ``nu_star_crit, Lambda_A^2 + Lambda_B^2, Lambda_alpha^2,
    Gamma_well`` per surface.  The offsets enter squared because their plain
    sum is unbounded below."""
    v = values
    return jnp.stack([v[:, I_OUT["nu_star_crit"]],
                      v[:, I_OUT["Lambda_A"]] ** 2 + v[:, I_OUT["Lambda_B"]] ** 2,
                      v[:, I_OUT["Lambda_alpha"]] ** 2,
                      v[:, I_OUT["Gamma_well"]]], axis=1)


def alpha_reduce(weights: AlphaWeights = AlphaWeights(),
                 surface_weights: Sequence[float] | None = None) -> Callable:
    """``F = sum_s w_s [w1 nu*_crit + w2 (Lambda_A^2 + Lambda_B^2) + w2b Lambda_alpha^2
    + w3 Gamma_well]``, each term divided by its reference."""
    w = jnp.asarray([weights.nu_star_crit, weights.offset, weights.Lambda_alpha,
                     weights.Gamma_well]) / jnp.asarray(weights.reference)

    def reduce(values):
        ws = (jnp.ones(values.shape[0]) if surface_weights is None
              else jnp.asarray(surface_weights, dtype=float))
        return jnp.sum(ws[:, None] * alpha_terms(values) * w[None, :])
    return reduce


def boundary_geometry(boundary, nfp: int, ntor: int, n_theta: int = 64, n_phi: int = 64):
    """``(volume, aspect)`` of the boundary ``[rbc, zbs]`` in VMEC's definition.

    ``R = sum rbc cos(m theta - n nfp phi)``, ``Z = sum zbs sin(...)``; the
    cross-section area ``A(phi) = |oint R dZ|``, the volume ``|int dphi oint
    R^2/2 dZ|``, ``Aminor = sqrt(<A>/pi)``, ``Rmajor = V/(2 pi <A>)``.  The
    trapezoidal rule is exact for the band-limited integrands once ``n_theta >
    3 mpol`` and ``n_phi > 3 nfp ntor``.
    """
    b = jnp.asarray(boundary)
    mpol = b.shape[1]
    m = jnp.arange(mpol)[:, None]
    n = (jnp.arange(b.shape[2]) - ntor)[None, :]
    th = jnp.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    ph = jnp.linspace(0.0, 2.0 * np.pi / nfp, n_phi, endpoint=False)
    arg = m[None, None] * th[:, None, None, None] - n[None, None] * nfp * ph[None, :, None, None]
    R = jnp.sum(b[0] * jnp.cos(arg), axis=(2, 3))
    dZ = jnp.sum(b[1] * m * jnp.cos(arg), axis=(2, 3))
    area = jnp.abs(jnp.mean(R * dZ, axis=0) * 2.0 * np.pi)
    volume = jnp.abs(jnp.mean(0.5 * R**2 * dZ) * 4.0 * np.pi**2)
    mean_area = jnp.mean(area)
    aminor = jnp.sqrt(mean_area / np.pi)
    rmajor = volume / (2.0 * np.pi * mean_area)
    return volume, rmajor / aminor


def aspect_penalty(target: float, weight: float, nfp: int, ntor: int, upper_only: bool = True):
    """``weight * max(A - target, 0)^2`` (or two-sided) of the boundary aspect ratio."""
    def c(boundary, harm):
        a = boundary_geometry(boundary, nfp, ntor)[1]
        d = jnp.maximum(a - target, 0.0) if upper_only else a - target
        return weight * d**2
    return c


def volume_penalty(target: float, weight: float, nfp: int, ntor: int):
    """``weight * (V/target - 1)^2`` of the plasma volume."""
    def c(boundary, harm):
        return weight * (boundary_geometry(boundary, nfp, ntor)[0] / target - 1.0) ** 2
    return c


def iota_penalty(target, weight: float):
    """``weight * mean((iotaf - target)^2)`` on the full mesh; ``target`` is an
    array over the radial grid or a scalar."""
    def c(boundary, harm):
        return weight * jnp.mean((harm["iotaf"] - jnp.asarray(target)) ** 2)
    return c


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
        reduce: Scalar of the ``(S, 7)`` outputs; default :func:`offset_reduce`.
        nu_star: Collisionality of the default offset objective.
        squared: Default objective ``sum_s L_s^2`` instead of ``sum_s L_s``.
        constraints: Penalties ``c(boundary, harm)`` (JAX scalars) added to the
            objective, with ``harm`` the primary wout arrays of the equilibrium.

    ``boundary`` is an array of shape ``(2, mpol, 2 ntor + 1)`` with rows
    ``rbc`` and ``zbs`` in VMEC++'s layout.  The last evaluation is cached, so
    ``objective`` followed by ``grad`` at the same boundary costs one solve.
    """

    def __init__(self, indata, s_list, settings: sd.DerivativeSettings, mboz: int, nboz: int,
                 reduce: Callable | None = None, nu_star: float = 1e-3, squared: bool = False,
                 workdir=None, constraints: Sequence[Callable] = ()):
        if indata.lasym:
            raise ValueError("stellarator-symmetric equilibria only")
        self.indata = indata.model_copy(update={"always_fix_m1_gauge": True})
        self.s = np.atleast_1d(np.asarray(s_list, dtype=float))
        self.settings = settings
        self.mboz, self.nboz = int(mboz), int(nboz)
        self.reduce = offset_reduce(nu_star, squared=squared) if reduce is None else reduce
        self._dreduce = jax.jit(jax.grad(self.reduce))
        self.constraints = tuple(constraints)

        def penalty(boundary, harm):
            return sum((c(boundary, harm) for c in self.constraints), jnp.float64(0.0))
        self._penalty = jax.jit(jax.value_and_grad(penalty, argnums=(0, 1)))
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
        """``(S, 7)`` outputs with the structure found at this boundary.

        The solve keeps the VMEC++ pullback, so a following :meth:`grad` at the
        same boundary needs only the adjoint solve.
        """
        return self.evaluate(boundary, with_pullback=True).values

    def penalty(self, boundary) -> float:
        """Sum of the constraint penalties (0 without constraints)."""
        if not self.constraints:
            return 0.0
        ev = self.evaluate(boundary, with_pullback=True)
        return float(self._penalty(jnp.asarray(ev.boundary), ev.harm)[0])

    def objective(self, boundary) -> float:
        return float(self.reduce(jnp.asarray(self.values(boundary)))) + self.penalty(boundary)

    def grad(self, boundary) -> np.ndarray:
        """Exact gradient of :meth:`objective`, shape ``(2, mpol, 2 ntor + 1)``."""
        ev = self.evaluate(boundary, with_pullback=True)
        return self.vjp(boundary, np.asarray(self._dreduce(jnp.asarray(ev.values))),
                        with_constraints=True)

    def vjp(self, boundary, cotangent, with_constraints: bool = False) -> np.ndarray:
        """``sum_{s,o} cotangent[s, o] d values[s, o] / d boundary``, plus the
        gradient of the constraint penalties with ``with_constraints``."""
        ev = self.evaluate(boundary, with_pullback=True)
        t0 = time.perf_counter()
        g_harm = ev.chain.vjp(np.asarray(cotangent, dtype=float), ev.harm)
        g_direct = 0.0
        if with_constraints and self.constraints:
            _, (gb, gh) = self._penalty(jnp.asarray(ev.boundary), ev.harm)
            g_harm = {k: np.asarray(g_harm[k]) + np.asarray(gh[k]) for k in g_harm}
            g_direct = np.asarray(gb)
        t1 = time.perf_counter()
        (g,) = ev.pullback({k: jnp.asarray(g_harm[k], dtype=ev.harm[k].dtype) for k in ev.harm})
        t2 = time.perf_counter()
        ev.timing.update({"pyrabe_vjp": t1 - t0, "vmecpp_adjoint": t2 - t1})
        return np.asarray(g) + g_direct

    def frozen_values(self, boundary, reference: Evaluation) -> np.ndarray:
        """Outputs at ``boundary`` with the discrete structure of ``reference``.

        This is the function whose derivative :meth:`grad` returns, so a finite
        difference of it checks the gradient without discrete switches.
        """
        harm = {k: jnp.asarray(np.asarray(v)) for k, v in self._harm_of_boundary(
            np.asarray(boundary, float)).items()}
        self.n_solves += 1
        return reference.chain.values(harm)
