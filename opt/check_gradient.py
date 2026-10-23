"""Finite-difference check of the boundary gradient over the full chain.

Every FD point re-solves VMEC++ (gauge pinned) and evaluates rabe's outputs by
the pyrabe primal, with the discrete structure frozen at the base boundary
(``frozen``) and, for the objective, also recomputed (``full``).  The gradient
side is VMEC++'s implicit adjoint composed with the JAX VJP of pyrabe.

usage: check_gradient.py INPUT NS S1,S2,.. M_POL N_TOR MBOZ NBOZ LINES NGAUSS [SEED [OUT.json]]
"""

from __future__ import annotations

import json
import sys
import time

import numpy as np
import vmecpp
from pyrabe import spectrum_derivative as sd

from rabe_boundary import BoundaryObjective


def load_input(path: str, ns: int, ftol: float = 1e-20, niter: int = 50000):
    inp = vmecpp.VmecInput.from_file(path)
    return inp.model_copy(update={"ns_array": np.array([ns]), "ftol_array": np.array([ftol]),
                                  "niter_array": np.array([niter])})


def smooth_direction(x0: np.ndarray, ntor: int, seed: int) -> np.ndarray:
    """Random direction scaled like the boundary spectrum, unused entries zero."""
    rng = np.random.default_rng(seed)
    _, mpol, nn = x0.shape
    m = np.arange(mpol)[:, None] * np.ones((1, nn))
    n = (np.arange(nn) - ntor)[None, :] * np.ones((mpol, 1))
    scale = np.abs(x0) + 1e-3 * np.abs(x0[0, 0, ntor]) * np.exp(-(m + np.abs(n)))
    d = rng.standard_normal(x0.shape) * scale
    d[:, (m == 0) & (n < 0)] = 0.0
    d[1, (m == 0) & (n == 0)] = 0.0
    return d / np.max(np.abs(d)) * np.max(np.abs(x0[:, 1:, :]))


def main(argv):
    path, ns, s_txt, M, N, mboz, nboz, lines, ng = argv[:9]
    seed = int(argv[9]) if len(argv) > 9 else 20261023
    s_list = [float(v) for v in s_txt.split(",")]
    st = sd.DerivativeSettings(M_pol=float(M), N_tor=float(N), max_n_fieldlines=int(lines),
                               n_gauss=int(ng))
    inp = load_input(path, int(ns))
    obj = BoundaryObjective(inp, s_list, st, int(mboz), int(nboz), nu_star=1e-3)
    x0 = obj.x0
    t = time.perf_counter()
    J0 = obj.objective(x0)
    g = obj.grad(x0)
    t_grad = time.perf_counter() - t
    ev = obj._last
    timing = dict(ev.timing, total_value_and_grad=t_grad)
    S = len(s_list)
    rows = {"objective": g}
    for o in (0, 1):
        for i in range(S):
            cot = np.zeros((S, 5))
            cot[i, o] = 1.0
            rows[f"{sd.OUTPUT_NAMES[o]}@{s_list[i]}"] = obj.vjp(x0, cot)
    d = smooth_direction(x0, int(inp.ntor), seed)
    exact = {k: float(np.sum(v * d)) for k, v in rows.items()}
    hs = [1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6]
    fd_frozen = {k: [] for k in rows}
    fd_full = []
    for h in hs:
        # A large step can move a line maximum out of reach of the frozen
        # structure (pyrabe locate raises); that step is recorded as NaN.
        try:
            vp = obj.frozen_values(x0 + h * d, ev)
            vm = obj.frozen_values(x0 - h * d, ev)
        except RuntimeError as err:
            print(f"h={h:.0e} frozen: {err}")
            vp = vm = np.full((S, 5), np.nan)
        dv = (vp - vm) / (2 * h)
        fd_frozen["objective"].append(float((obj.reduce(vp) - obj.reduce(vm)) / (2 * h)))
        for o in (0, 1):
            for i in range(S):
                fd_frozen[f"{sd.OUTPUT_NAMES[o]}@{s_list[i]}"].append(float(dv[i, o]))
        try:
            fd_full.append((obj.objective(x0 + h * d) - obj.objective(x0 - h * d)) / (2 * h))
        except RuntimeError as err:
            print(f"h={h:.0e} full: {err}")
            fd_full.append(float("nan"))
    out = {"input": path, "ns": int(ns), "s": s_list, "J0": J0, "values": ev.values.tolist(),
           "timing": timing, "h": hs, "exact": exact, "fd_frozen": fd_frozen, "fd_full": fd_full}
    print(f"ns={ns} s={s_list} J0={J0:.6e} timing={ {k: round(v, 2) for k, v in timing.items()} }")
    for k, e in exact.items():
        rel = [abs(f - e) / max(abs(e), 1e-300) for f in fd_frozen[k]]
        print(f"{k:22s} exact {e:+.6e}  frozen-FD rel err " + " ".join(f"{r:.1e}" for r in rel))
    rel = [abs(f - exact["objective"]) / abs(exact["objective"]) for f in fd_full]
    print(f"{'objective (full)':22s} exact {exact['objective']:+.6e}  full-FD rel err   "
          + " ".join(f"{r:.1e}" for r in rel))
    return out


if __name__ == "__main__":
    res = main(sys.argv[1:])
    if len(sys.argv) > 11:
        with open(sys.argv[11], "w") as f:
            json.dump(res, f, indent=1)
