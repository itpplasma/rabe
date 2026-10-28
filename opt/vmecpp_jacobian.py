"""Warm-started VMEC++ solves and batched implicit adjoints for full Jacobians.

``vmecpp.run`` with a traced boundary (VMEC++ 0.9.x) solves cold from the
multigrid schedule and, in its VJP, assembles and LU-factors the transposed
force Hessian anew for every cotangent; a second cotangent even re-solves the
equilibrium, because the forward state is popped from a one-shot cache.  A
least-squares optimizer needs the Jacobian of ``K`` residuals, which on that
path costs ``K`` solves and ``K`` factorizations.

:func:`traced_run` runs the same traced ``vmecpp.run`` with
:class:`_WarmRunSolver` in place of ``vmecpp.autodiff._RunSolver``:

- the forward solve hot-restarts from the last converged equilibrium of the
  same input (``restart_from`` of the untraced ``vmecpp.run``, final ``ns``
  only); if that fails or misses ``ftol`` it falls back to the cold multigrid
  schedule, so the converged state is that of the input's tolerance either way;
- the VJP keeps the solved model and the factorized block-tridiagonal ``H^T``
  of one state and applies them to every cotangent, so ``jax.vmap`` over the
  pullback of ``K`` cotangents costs one assembly, one factorization and ``K``
  block back-substitutions.

The adjoint arithmetic per cotangent is that of ``vmecpp.autodiff._implicit_vjp``
(same deflation, same residual check); only the factorization is shared.  This
uses private functions of vmecpp 0.9.1 and checks the version on import.
"""

from __future__ import annotations

import collections
import contextlib
import time

import jax
import numpy as np
import vmecpp
from vmecpp import _vmecpp
from vmecpp import autodiff as ad

SUPPORTED = ("0.9.",)
_version = getattr(vmecpp, "__version__", None)
if _version is None:
    from importlib.metadata import version as _pkg_version
    _version = _pkg_version("vmecpp")
if not _version.startswith(SUPPORTED):
    raise ImportError(f"vmecpp_jacobian relies on private vmecpp {SUPPORTED} internals, "
                      f"found {_version}")


class SolveStats:
    """Counters and the last solve's diagnostics, shared by all solvers."""

    def __init__(self):
        self.reset()

    def reset(self):
        self.forward = 0
        self.warm = 0
        self.warm_fallback = 0
        self.factorizations = 0
        self.adjoint_rows = 0
        self.t_forward = 0.0
        self.t_factor = 0.0
        self.t_rows = 0.0
        self.last: dict = {}

    def as_dict(self):
        return {k: v for k, v in self.__dict__.items()}


STATS = SolveStats()

# The last converged (cpp indata, cpp wout) per warm-start key, for hot restarts.
_WARM: dict[object, tuple] = {}
# Solved models and their factorized adjoint operators, keyed by resolution and boundary.
_ADJOINT: collections.OrderedDict[tuple, dict] = collections.OrderedDict()
_ADJOINT_SIZE = 2


def _converged(wout, ftol: float) -> bool:
    return (int(wout.ier_flag) == 0
            and max(float(wout.fsqr), float(wout.fsqz), float(wout.fsql)) <= ftol)


class _WarmRunSolver(ad._RunSolver):
    """``_RunSolver`` with a hot-restarted forward solve and a shared adjoint factorization."""

    warm_key: object = None  # set by traced_run for the duration of one call

    def _run(self, boundary: np.ndarray, profiles: dict):
        indata = ad._make_indata(self.vmec_input._to_cpp_vmecindata(), boundary, profiles)
        indata.always_fix_m1_gauge = True
        ftol = float(np.asarray(indata.ftol_array)[-1])
        key = type(self).warm_key
        t0 = time.perf_counter()
        output, warm = None, False
        if key is not None and key in _WARM:
            prev_indata, prev_wout = _WARM[key]
            single = indata.copy()
            single.ns_array = np.asarray(indata.ns_array)[-1:]
            single.ftol_array = np.asarray(indata.ftol_array)[-1:]
            single.niter_array = np.asarray(indata.niter_array)[-1:]
            try:
                out = _vmecpp.run(single, initial_state=_vmecpp.HotRestartState(
                    wout=prev_wout, indata=prev_indata), max_threads=self.max_threads,
                    verbose=_vmecpp.OutputMode(0))
                if _converged(out.wout, ftol):
                    output, warm, indata = out, True, single
                else:
                    STATS.warm_fallback += 1
            except Exception:  # noqa: BLE001 - fall back to the cold schedule below
                STATS.warm_fallback += 1
        if output is None:
            output = _vmecpp.run(indata, max_threads=self.max_threads,
                                 verbose=_vmecpp.OutputMode(self.verbose))
        STATS.forward += 1
        STATS.warm += int(warm)
        STATS.t_forward += time.perf_counter() - t0
        w = output.wout
        STATS.last = dict(warm=warm, fsqr=float(w.fsqr), fsqz=float(w.fsqz), fsql=float(w.fsql),
                          itfsq=int(w.itfsq), ier_flag=int(w.ier_flag), ftol=ftol,
                          converged=_converged(w, ftol), t=time.perf_counter() - t0)
        if key is not None and STATS.last["converged"]:
            _WARM[key] = (indata, w)
        initial_state = _vmecpp.HotRestartState(wout=w, indata=indata)
        model = _vmecpp.VmecModel.create(indata, self.ns, initial_state=initial_state)
        flip = -1.0 if model.have_to_flip_theta else 1.0
        return output, np.array(model.get_state()), flip

    def _backward_callback(self, boundary: np.ndarray, geometry_bar: np.ndarray, profiles=None):
        if profiles:
            return super()._backward_callback(boundary, geometry_bar, profiles)
        key = (self.ns, self.vmec_input.mpol, self.vmec_input.ntor,
               repr(np.asarray(self.vmec_input.ftol_array).tolist()),
               np.ascontiguousarray(boundary, dtype=np.float64).tobytes())
        entry = _ADJOINT.get(key)
        if entry is None:
            t0 = time.perf_counter()
            model = self._solved_model(boundary, {})
            entry = _factorize(model)
            _ADJOINT[key] = entry
            while len(_ADJOINT) > _ADJOINT_SIZE:
                _ADJOINT.popitem(last=False)
            STATS.factorizations += 1
            STATS.t_factor += time.perf_counter() - t0
        t0 = time.perf_counter()
        bar = _apply_adjoint(entry, np.asarray(geometry_bar, dtype=np.float64))
        STATS.adjoint_rows += 1
        STATS.t_rows += time.perf_counter() - t0
        return bar, np.zeros((3, self.ns - 1))


def _factorize(model) -> dict:
    """The parts of ``_implicit_vjp`` that depend only on the solved state."""
    state = np.asarray(model.get_state(), dtype=np.float64)
    interior, boundary = ad._interior_and_boundary(model)
    gauge = ad._gauge_entries(model)
    interior = np.setdiff1d(interior, gauge)
    prescribed = np.concatenate([boundary, gauge])
    model.set_state(np.ascontiguousarray(state))
    model.evaluate(2, 2, True)
    interior = ad._structural_nullfree_interior(model, interior)
    matrix = ad._assemble_block_tridiagonal(model, interior)
    matrix.factorize()
    return dict(model=model, interior=interior, prescribed=prescribed, gauge=gauge,
                matrix=matrix, state_size=state.size)


def _apply_adjoint(entry: dict, geometry_bar: np.ndarray) -> np.ndarray:
    """``_implicit_vjp(model, geometry_bar)[0]`` with the factorization of ``entry``."""
    model, interior, prescribed = entry["model"], entry["interior"], entry["prescribed"]
    coefficient_bar = np.asarray(geometry_bar[2 * model.ns:], dtype=np.float64)
    poloidal_flux_bar = np.asarray(geometry_bar[model.ns:2 * model.ns], dtype=np.float64)
    state_bar = np.asarray(model.geometry_state_vjp(coefficient_bar, poloidal_flux_bar),
                           dtype=np.float64)
    embedded = np.zeros(entry["state_size"])
    embedded[interior] = entry["matrix"].solve(state_bar[interior])
    transposed = np.asarray(model.exact_hessian_vector_product_transpose(
        np.ascontiguousarray(embedded)), dtype=np.float64)
    residual = float(np.linalg.norm(transposed[interior] - state_bar[interior]))
    scale = max(float(np.linalg.norm(state_bar[interior])), 1.0e-300)
    if residual > 1.0e-6 * scale:
        raise RuntimeError(f"VMEC++ implicit adjoint solve failed: relative residual "
                           f"{residual / scale:.3e}")
    full = np.zeros(entry["state_size"])
    full[prescribed] = state_bar[prescribed] - transposed[prescribed]
    full = ad._pull_gauge_back_to_boundary(model, full, entry["gauge"])
    return ad._boundary_from_state_vjp(model, full)


@contextlib.contextmanager
def _patched(warm_key):
    original = ad._RunSolver
    ad._RunSolver = _WarmRunSolver
    _WarmRunSolver.warm_key = warm_key
    try:
        yield
    finally:
        ad._RunSolver = original
        _WarmRunSolver.warm_key = None


def traced_run(indata, boundary, *, max_threads=None, warm_key=None):
    """``vmecpp.run`` of ``indata`` with the traced ``boundary = [rbc, zbs]``.

    Use inside ``jax.vjp``.  ``warm_key`` (any hashable, for example the
    optimization run) enables hot restarts from the last converged solve with
    the same key; ``None`` solves cold.
    """
    with _patched(warm_key):
        return vmecpp.run(indata.model_copy(update={"rbc": boundary[0], "zbs": boundary[1]}),
                          max_threads=max_threads, verbose=False)


def forget_warm(warm_key=None):
    """Drop the warm-start state of ``warm_key`` (all keys with ``None``)."""
    if warm_key is None:
        _WARM.clear()
    else:
        _WARM.pop(warm_key, None)


def batched_pullback(pullback, cotangents):
    """``[pullback(c) for c in cotangents]`` as one ``jax.vmap`` call.

    ``cotangents`` is a pytree like the primal output with a leading batch axis.
    With :func:`traced_run` the VMEC++ adjoint is factorized once for the batch.
    """
    (g,) = jax.vmap(pullback)(cotangents)
    return np.asarray(g)
