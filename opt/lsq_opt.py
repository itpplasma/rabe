"""Least-squares boundary optimization with the exact Jacobian of rabe's criteria.

The objective is written as ``F(x) = |r(x)|^2`` with a residual vector ``r``:

- per surface and objective term (:data:`OBJECTIVES`), each term normalized so
  that its sum over the surfaces is 1 at the start boundary:
  ``Lambda_A``, ``Lambda_B``, ``Lambda_alpha`` enter linearly (their squares are
  the terms), ``nu_star_crit``, ``Gamma_well``, ``var_dJ``, ``Gamma_s`` through
  ``sqrt(v + floor)`` with ``floor = 1e-8`` of the start value (their sums are
  the terms, up to that constant, and the square root stays differentiable);
- constraints ``c(boundary, harm)``: equalities (aspect ratio and volume
  relative to the start) and inequalities ``c <= 0`` (``iotaf`` inside the start
  range widened by ``iota_margin``, the Boozer mirror ratio of every objective
  surface inside ``1 +- mirror_band`` of its start value), entering either as
  quadratic penalties ``sqrt(w) c`` / ``sqrt(w) max(c, 0)`` or, for the
  augmented-Lagrangian solver, as ``sqrt(mu) (c + lambda/mu)`` and
  ``sqrt(mu) max(c + lambda/mu, 0)``.

``F`` with penalties equals ``alpha_optimize``'s objective for E1 and E2 (up to
the floor constant and the mirror penalty normalized per surface).

Jacobian of one evaluation (:meth:`Problem.evaluate`): one VMEC++ solve under
``jax.vjp`` (hot-restarted from the previous solve, :mod:`vmecpp_jacobian`);
one rabe structure and values; ``dR/d(values, boundary, harm)`` by ``jax.jacrev``
of the (jitted) residual map; per surface one VJP of the pyrabe chain per used
output; and one batched VMEC++ implicit adjoint (one factorization, one
back-substitution per residual that depends on the equilibrium).  The gradient
mode needs a single adjoint row.  rabe's discrete structure is frozen at every
evaluated boundary, as in :mod:`rabe_boundary`.
"""

from __future__ import annotations

import os
import time
from collections.abc import Sequence
from dataclasses import dataclass, field

import jax
import jax.numpy as jnp
import numpy as np
import vmecpp
from pyrabe import boozer_libneo as bl
from pyrabe import spectrum_derivative as sd
from pyrabe.wout_derivative import WoutChain, spectrum_from_harmonics

import vmecpp_jacobian as vj
from rabe_boundary import boundary_geometry, harmonics_of, mirror_ratio

jax.config.update("jax_enable_x64", True)

# term -> (outputs, form); "lin": the outputs are residuals, "sqrt": sqrt(v + floor)
TERMS = {
    "nu_star_crit": (("nu_star_crit",), "sqrt"),
    "offset": (("Lambda_A", "Lambda_B"), "lin"),
    "Lambda_alpha": (("Lambda_alpha",), "lin"),
    "Gamma_well": (("Gamma_well",), "sqrt"),
    "var_dJ": (("var_dJ",), "sqrt"),
    "Gamma_s": (("Gamma_s",), "sqrt"),
}
OBJECTIVES = {
    "E1": ("nu_star_crit", "offset"),
    "E2": ("nu_star_crit", "Lambda_alpha", "Gamma_well"),
    "E3": ("nu_star_crit", "var_dJ", "Gamma_s"),
}
FLOOR = 1e-8
BOOZER_ROUTE = "pyrabe.boozer_libneo (libneo hybrid transform from VMEC++ primaries rmnc, zmns, lmns, iotaf)"


def term_values(values, objective: str):
    """``(S, T)`` term values per surface (squares for "lin", plain for "sqrt")."""
    idx = {n: i for i, n in enumerate(sd.OUTPUT_NAMES)}
    cols = []
    for t in OBJECTIVES[objective]:
        outs, form = TERMS[t]
        missing = [o for o in outs if o not in idx]
        if missing:
            raise KeyError(f"objective {objective} needs pyrabe outputs {missing}; "
                           f"this pyrabe has {sd.OUTPUT_NAMES}")
        v = jnp.stack([values[:, idx[o]] for o in outs], axis=1)
        cols.append(jnp.sum(v**2, axis=1) if form == "lin" else v[:, 0])
    return jnp.stack(cols, axis=1)


@dataclass
class Reference:
    """Start-point normalization and constraint targets (one per ``ns``)."""

    term_sums: np.ndarray       # (T,)
    aspect: float
    volume: float
    iota_min: float
    iota_max: float
    mirror: np.ndarray          # (S,)

    def to_json(self):
        return dict(term_sums=self.term_sums.tolist(), aspect=self.aspect, volume=self.volume,
                    iota_min=self.iota_min, iota_max=self.iota_max, mirror=self.mirror.tolist())

    @classmethod
    def from_json(cls, d):
        return cls(np.asarray(d["term_sums"]), d["aspect"], d["volume"], d["iota_min"],
                   d["iota_max"], np.asarray(d["mirror"]))


@dataclass
class Weights:
    w_aspect: float = 1e3
    w_volume: float = 1e3
    w_iota: float = 1e3
    w_mirror: float = 1e3
    iota_margin: float = 0.02
    mirror_band: float = 0.1


@dataclass
class Evaluation:
    x: np.ndarray
    boundary: np.ndarray
    ok: bool
    r: np.ndarray | None = None
    J: np.ndarray | None = None          # (K, n) in x
    g: np.ndarray | None = None          # gradient of |r|^2 in x
    values: np.ndarray | None = None
    terms: np.ndarray | None = None
    c_eq: np.ndarray | None = None
    c_in: np.ndarray | None = None
    info: dict = field(default_factory=dict)

    @property
    def F(self):
        return float(self.r @ self.r) if self.ok else np.inf

    @property
    def objective(self):
        return float(np.sum(self.terms)) if self.terms is not None else np.inf

    def violation(self):
        if not self.ok:
            return np.inf
        return float(max(np.max(np.abs(self.c_eq), initial=0.0), np.max(self.c_in, initial=0.0)))


class Problem:
    """Residuals and Jacobian for one VMEC++ resolution and one set of free coefficients.

    Args:
        indata: VMEC++ input with the resolution of this stage (``ns_array``, ``ftol_array``).
        s_list, settings, mboz, nboz: rabe surfaces and settings.
        objective: ``"E1"``, ``"E2"`` or ``"E3"``.
        mask: Free entries of ``[rbc, zbs]``; ``x = (b - b0)[mask] / scale``.
        b0, scale: Reference boundary and the length of ``x = 1``.
        weights: Penalty weights and constraint bands.
        reference: Start normalization (computed at ``b0`` when ``None``).
        workdir: Holds the base wout of this resolution (libneo plans).
        warm_key: Hot-restart key of :mod:`vmecpp_jacobian` (``None``: cold solves).
    """

    def __init__(self, indata, s_list, settings, mboz, nboz, objective, mask, b0, scale,
                 weights: Weights, workdir, reference: Reference | None = None,
                 warm_key=None, max_threads=None):
        self.inp = indata.model_copy(update={"always_fix_m1_gauge": True,
                                             "return_outputs_even_if_not_converged": True})
        self.s = np.atleast_1d(np.asarray(s_list, float))
        self.settings, self.mboz, self.nboz = settings, int(mboz), int(nboz)
        self.objective = objective
        self.mask, self.b0, self.scale = mask, np.asarray(b0, float), float(scale)
        self.weights = weights
        self.nfp, self.ntor = int(self.inp.nfp), int(self.inp.ntor)
        self.ns = int(np.asarray(self.inp.ns_array)[-1])
        self.ftol = float(np.asarray(self.inp.ftol_array)[-1])
        self.warm_key, self.max_threads = warm_key, max_threads
        os.makedirs(workdir, exist_ok=True)
        self.base_wout = os.path.join(workdir, f"wout_base_ns{self.ns}.nc")
        if not os.path.exists(self.base_wout):
            vmecpp.run(self.inp, max_threads=max_threads, verbose=False).wout.save(self.base_wout)
        self.plans = [bl.cached_plan(self.base_wout, float(s)) for s in self.s]
        self.n_eval = 0
        self.mu, self.lam_eq, self.lam_in = None, None, None   # augmented Lagrangian state
        self.reference = reference
        self._jit = {}
        if self.reference is None:
            ev = self._forward(self.b0)
            if not ev.ok:
                raise RuntimeError(f"start boundary failed: {ev.info.get('error')}")
            self.reference = self._reference(ev)

    # -- maps ---------------------------------------------------------------------------
    @property
    def n(self):
        return int(self.mask.sum())

    def boundary(self, x):
        b = self.b0.copy()
        b[self.mask] += self.scale * np.asarray(x, float)
        return b

    def x_of(self, boundary):
        return (np.asarray(boundary, float) - self.b0)[self.mask] / self.scale

    def mirror(self, harm):
        out = []
        for plan in self.plans:
            p, xm, xn = spectrum_from_harmonics(plan, harm, self.mboz, self.nboz)
            out.append(mirror_ratio(p.bmnc, xm, xn))
        return jnp.stack(out)

    def constraints(self, boundary, harm):
        """``(c_eq, c_in)``: relative aspect and volume; iota and mirror bands (``<= 0``)."""
        ref, w = self.reference, self.weights
        vol, asp = boundary_geometry(boundary, self.nfp, self.ntor)
        c_eq = jnp.stack([asp / ref.aspect - 1.0, vol / ref.volume - 1.0])
        iota = harm["iotaf"]
        d = self.mirror(harm)
        d0 = jnp.asarray(ref.mirror)
        c_in = jnp.concatenate([iota - (ref.iota_max + w.iota_margin),
                                (ref.iota_min - w.iota_margin) - iota,
                                (d - d0 * (1 + w.mirror_band)) / d0,
                                (d0 * (1 - w.mirror_band) - d) / d0])
        return c_eq, c_in

    def constraint_weights(self):
        w, S = self.weights, self.s.size
        n_iota = self.ns
        w_eq = np.array([w.w_aspect, w.w_volume])
        w_in = np.concatenate([np.full(2 * n_iota, w.w_iota / n_iota), np.full(2 * S, w.w_mirror / S)])
        return w_eq, w_in

    def residual_fn(self, values, boundary, harm, params):
        """Residual vector; ``params`` holds the normalization and the multipliers."""
        terms = term_values(values, self.objective)                       # (S, T)
        cols = []
        for j, t in enumerate(OBJECTIVES[self.objective]):
            outs, form = TERMS[t]
            scale = jnp.sqrt(1.0 / params["term_sums"][j])
            if form == "lin":
                idx = [sd.OUTPUT_NAMES.index(o) for o in outs]
                cols.append((values[:, idx] * scale).ravel())
            else:
                cols.append(jnp.sqrt(terms[:, j] + params["floor"][j]) * scale)
        c_eq, c_in = self.constraints(boundary, harm)
        if params["mode"] == "al":
            mu = params["mu"]
            pen = [jnp.sqrt(mu) * (c_eq + params["lam_eq"] / mu),
                   jnp.sqrt(mu) * jnp.maximum(c_in + params["lam_in"] / mu, 0.0)]
        else:
            pen = [jnp.sqrt(params["w_eq"]) * c_eq,
                   jnp.sqrt(params["w_in"]) * jnp.maximum(c_in, 0.0)]
        return jnp.concatenate(cols + pen)

    def params(self):
        ref = self.reference
        w_eq, w_in = self.constraint_weights()
        p = dict(term_sums=jnp.asarray(ref.term_sums),
                 floor=jnp.asarray(FLOOR * ref.term_sums / self.s.size),
                 w_eq=jnp.asarray(w_eq), w_in=jnp.asarray(w_in), mode="pen")
        if self.mu is not None:
            p.update(mode="al", mu=self.mu, lam_eq=jnp.asarray(self.lam_eq),
                     lam_in=jnp.asarray(self.lam_in))
        return p

    def _residual_jit(self, mode):
        key = ("r", mode)
        if key not in self._jit:
            def f(values, boundary, harm, params):
                return self.residual_fn(values, boundary, harm, dict(params, mode=mode))
            self._jit[key] = jax.jit(f)
            self._jit[("J", mode)] = jax.jit(jax.jacrev(f, argnums=(0, 1, 2)))
        return self._jit[key], self._jit[("J", mode)]

    # -- evaluation ---------------------------------------------------------------------
    def _harm(self, boundary):
        return harmonics_of(vj.traced_run(self.inp, boundary, max_threads=self.max_threads,
                                          warm_key=self.warm_key).wout)

    def _forward(self, boundary, x=None):
        t0 = time.perf_counter()
        ev = Evaluation(x=None if x is None else np.asarray(x, float), boundary=boundary, ok=False)
        try:
            harm, pullback = jax.vjp(self._harm, jnp.asarray(boundary))
            solve = dict(vj.STATS.last)
            if not solve.get("converged", False):
                raise RuntimeError(f"VMEC++ not converged: {solve}")
            harm = {k: jnp.asarray(np.asarray(v)) for k, v in harm.items()}
            t1 = time.perf_counter()
            chain = WoutChain(self.base_wout, self.s, self.settings, self.mboz, self.nboz, harm=harm)
            values = chain.values(harm)
            if not np.all(np.isfinite(values)):
                raise FloatingPointError("non-finite rabe values")
            ev.values, ev.ok = values, True
            ev.info.update(solve=solve, t_solve=t1 - t0, t_values=time.perf_counter() - t1)
            ev._harm, ev._pullback, ev._chain = harm, pullback, chain
        except Exception as e:  # noqa: BLE001 - a failed evaluation is a rejected step
            ev.info["error"] = f"{type(e).__name__}: {e}"
        return ev

    def _reference(self, ev):
        terms = np.asarray(term_values(jnp.asarray(ev.values), self.objective))
        iota = np.asarray(ev._harm["iotaf"])
        vol, asp = (float(v) for v in boundary_geometry(self.b0, self.nfp, self.ntor))
        sums = terms.sum(axis=0)
        return Reference(np.where(sums > 0, sums, 1.0), asp, vol, float(iota.min()),
                         float(iota.max()), np.asarray(self.mirror(ev._harm)))

    def evaluate(self, x, mode: str = "r") -> Evaluation:
        """``mode``: ``"r"`` residuals, ``"J"`` plus the Jacobian, ``"g"`` plus ``grad |r|^2``."""
        x = np.asarray(x, float)
        b = self.boundary(x)
        last = getattr(self, "_last_fwd", None)
        if last is not None and np.array_equal(last.boundary, b):
            ev = Evaluation(x=x, boundary=b, ok=last.ok, values=last.values, info=dict(last.info))
            ev.info["reused_forward"] = True
            if last.ok:
                ev._harm, ev._pullback, ev._chain = last._harm, last._pullback, last._chain
        else:
            ev = self._forward(b, x)
            self.n_eval += 1
            self._last_fwd = ev if ev.ok else None
        if not ev.ok:
            return ev
        params = self.params()
        mode_p = params.pop("mode")
        f, jac = self._residual_jit(mode_p)
        harm, values = ev._harm, jnp.asarray(ev.values)
        ev.r = np.asarray(f(values, jnp.asarray(b), harm, params))
        ev.terms = np.asarray(term_values(values, self.objective))
        ev.terms = ev.terms / self.reference.term_sums[None, :]
        c_eq, c_in = self.constraints(jnp.asarray(b), harm)
        ev.c_eq, ev.c_in = np.asarray(c_eq), np.asarray(c_in)
        if not np.all(np.isfinite(ev.r)):
            ev.ok, ev.info["error"] = False, "non-finite residual"
            return ev
        if mode in ("J", "g"):
            t0 = time.perf_counter()
            dv, db, dh = jac(values, jnp.asarray(b), harm, params)
            dv, db = np.asarray(dv), np.asarray(db)
            dh = {k: np.asarray(v) for k, v in dh.items()}
            if mode == "g":   # contract with 2 r first: one adjoint row
                w = 2.0 * ev.r
                dv, db = np.tensordot(w, dv, 1)[None], np.tensordot(w, db, 1)[None]
                dh = {k: np.tensordot(w, v, 1)[None] for k, v in dh.items()}
            K = dv.shape[0]
            # harm cotangents of every row: direct part plus the pyrabe chain
            cot = {k: dh[k].copy() for k in dh}
            for i, (_, _, fs) in enumerate(ev._chain._prepare(harm)):
                rows = dv[:, i, :]                                   # (K, N_OUT)
                used = np.flatnonzero(np.any(rows != 0.0, axis=0))
                if used.size == 0:
                    continue
                _, pb = jax.vjp(fs, harm)
                if mode == "g":
                    (gi,) = pb(jnp.asarray(rows[0]))
                    for k in cot:
                        cot[k][0] += np.asarray(gi[k])
                else:
                    eye = jnp.eye(len(sd.OUTPUT_NAMES))[used]
                    (G,) = jax.vmap(pb)(eye)                          # (U, *harm)
                    for k in cot:
                        cot[k] += np.tensordot(rows[:, used], np.asarray(G[k]), 1)
            t1 = time.perf_counter()
            active = np.flatnonzero([any(np.any(cot[k][r] != 0.0) for k in cot) for r in range(K)])
            Jb = db.reshape(K, *b.shape).copy()
            if active.size:
                batch = {k: jnp.asarray(cot[k][active], dtype=harm[k].dtype) for k in cot}
                Jb[active] += vj.batched_pullback(ev._pullback, batch)
            t2 = time.perf_counter()
            Jx = Jb[:, self.mask] * self.scale
            if not np.all(np.isfinite(Jx)):
                ev.ok, ev.info["error"] = False, "non-finite Jacobian"
                return ev
            if mode == "g":
                ev.g = Jx[0]
            else:
                ev.J = Jx
                ev.g = 2.0 * Jx.T @ ev.r
            ev.info.update(t_pyrabe=t1 - t0, t_adjoint=t2 - t1, adjoint_rows=int(active.size))
        if ev is not self._last_fwd:
            for k in ("_harm", "_pullback", "_chain"):
                ev.__dict__.pop(k, None)
        return ev


def summary(ev: Evaluation) -> dict:
    out = dict(ok=ev.ok, info=ev.info)
    if ev.ok:
        out.update(F=ev.F, objective=ev.objective, violation=ev.violation(),
                   terms=ev.terms.tolist(), values=ev.values.tolist(),
                   c_eq=ev.c_eq.tolist(), c_in_max=float(np.max(ev.c_in)))
    return out


# -- solvers ----------------------------------------------------------------------------

def levenberg_marquardt(problem: Problem, x0, *, maxiter=50, lam0=1e-2, max_step=5.0, bound=50.0,
                        ftol=1e-10, xtol=1e-8, on_iter=None, on_eval=None, state=None):
    """Levenberg-Marquardt in ``x`` with an explicit step cap ``|dx|_2 <= max_step``.

    ``(J^T J + lam I) dx = -J^T r``; a trial is accepted when the gain ratio
    ``rho = (F - F_trial) / (F - |r + J dx|^2) > 1e-3``; ``lam`` follows
    Nielsen's update.  Trials need residuals only, the Jacobian is evaluated at
    accepted points.  ``state`` restores ``(x, lam, nu, k)`` from a checkpoint.
    """
    st = dict(x=np.asarray(x0, float), lam=lam0, nu=2.0, k=0, nfev=0, njev=0)
    if state:
        st.update(state)
    ev = problem.evaluate(st["x"], "J")
    st["nfev"] += 1
    st["njev"] += 1
    if on_eval:
        on_eval(ev, "J")
    if not ev.ok:
        raise RuntimeError(f"LM start failed: {ev.info.get('error')}")
    if on_iter:
        on_iter(ev, st)
    while st["k"] < maxiter:
        J, r = ev.J, ev.r
        A, g = J.T @ J, J.T @ r
        accepted = False
        for _ in range(12):
            dx = np.linalg.solve(A + st["lam"] * np.eye(A.shape[0]), -g)
            nrm = np.linalg.norm(dx)
            if nrm > max_step:
                dx *= max_step / nrm
            xt = np.clip(st["x"] + dx, -bound, bound)
            dx = xt - st["x"]
            pred = ev.F - float(np.sum((r + J @ dx) ** 2))
            et = problem.evaluate(xt, "r")
            st["nfev"] += 1
            if on_eval:
                on_eval(et, "r")
            rho = (ev.F - et.F) / pred if (et.ok and pred > 0) else -1.0
            if rho > 1e-3:
                st["lam"] *= max(1.0 / 3.0, 1.0 - (2.0 * rho - 1.0) ** 3)
                st["nu"] = 2.0
                accepted = True
                break
            st["lam"] *= st["nu"]
            st["nu"] *= 2.0
        if not accepted:
            return ev, st, "no acceptable step"
        ev_new = problem.evaluate(xt, "J")
        st["nfev"] += 1
        st["njev"] += 1
        if on_eval:
            on_eval(ev_new, "J")
        if not ev_new.ok:  # residual was fine but the Jacobian failed
            return ev, st, f"Jacobian failed: {ev_new.info.get('error')}"
        dF = ev.F - ev_new.F
        ev, st["x"], st["k"] = ev_new, xt, st["k"] + 1
        if on_iter:
            on_iter(ev, st)
        if dF <= ftol * max(ev.F, 1e-300) or np.linalg.norm(dx) <= xtol * (1 + np.linalg.norm(xt)):
            return ev, st, "converged"
    return ev, st, "maxiter"


def augmented_lagrangian(problem: Problem, x0, *, outer=6, inner=8, mu0=1e3, mu_growth=10.0,
                         tol=1e-3, on_iter=None, on_eval=None, state=None, **lm):
    """Equality and inequality constraints by an augmented Lagrangian with LM inner solves."""
    n_eq, n_in = 2, 2 * problem.ns + 2 * problem.s.size
    st = dict(outer=0, mu=mu0, lam_eq=np.zeros(n_eq), lam_in=np.zeros(n_in), x=np.asarray(x0, float),
              viol=np.inf, lm=None)
    if state:
        st.update({k: (np.asarray(v) if k.startswith("lam") else v) for k, v in state.items()})
    ev, msg = None, ""
    while st["outer"] < outer:
        problem.mu, problem.lam_eq, problem.lam_in = st["mu"], st["lam_eq"], st["lam_in"]
        def inner_iter(ev, s):
            if on_iter:
                on_iter(ev, {**{k: st[k] for k in ("outer", "mu", "lam_eq", "lam_in", "viol")},
                             "lm": {k: s[k] for k in ("lam", "nu", "k")}, "k": s["k"]})
        ev, lmst, msg = levenberg_marquardt(problem, st["x"], maxiter=inner, on_iter=inner_iter,
                                            on_eval=on_eval, state=st["lm"], **lm)
        st["x"] = lmst["x"]
        st["lm"] = None
        viol = ev.violation()
        st["lam_eq"] = st["lam_eq"] + st["mu"] * ev.c_eq
        st["lam_in"] = np.maximum(0.0, st["lam_in"] + st["mu"] * ev.c_in)
        if viol > 0.25 * st["viol"]:
            st["mu"] *= mu_growth
        st["viol"] = viol
        st["outer"] += 1
        if on_iter:
            on_iter(ev, dict({k: v for k, v in st.items() if k != "lm"}, al=True))
        if viol <= tol and msg == "converged":
            break
    problem.mu = None
    return ev, st, msg
