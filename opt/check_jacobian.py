"""Finite-difference V-curve of the least-squares Jacobian of :mod:`lsq_opt`.

Central differences of the residual vector along one random unit direction in
``x``, for each step ``h``, against ``J d``:

- ``refound``: each FD point re-solves VMEC++ and re-finds rabe's structure
  (the function the optimizer sees);
- ``frozen``: each FD point re-solves VMEC++ but keeps the structure of the
  base point (the function whose derivative ``J`` is).

Errors are reported per residual group (objective terms, equalities, iota and
mirror inequalities) relative to the largest ``|J d|`` of that group.  Also
checks the gradient mode against ``2 J^T r`` and times one Jacobian.

usage: check_jacobian.py INPUT NS OBJECTIVE M_POL N_TOR S1,S2 FREE_M OUT.json
       [--lines 60 --n-gauss 32 --mboz 9 --nboz 8 --ftol 1e-20 --x-scale 0.3 --threads 4]
"""

from __future__ import annotations

import argparse
import json
import os
import tempfile
import time

import jax.numpy as jnp
import numpy as np
import vmecpp

import lsq_opt as lo
from alpha_optimize import free_mask, settings_of


def main():
    ap = argparse.ArgumentParser()
    for a in ("input", "ns", "objective", "M_pol", "N_tor", "s", "free_m", "out"):
        ap.add_argument(a)
    ap.add_argument("--lines", type=int, default=60)
    ap.add_argument("--n-gauss", type=int, default=32)
    ap.add_argument("--mboz", type=int, default=9)
    ap.add_argument("--nboz", type=int, default=8)
    ap.add_argument("--ftol", type=float, default=1e-20)
    ap.add_argument("--niter", type=int, default=100000)
    ap.add_argument("--x-scale", type=float, default=0.3)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--hs", default="1,0.3,0.1,0.03,0.01,3e-3,1e-3")
    a = ap.parse_args()
    ns = int(a.ns)
    inp = vmecpp.VmecInput.from_file(a.input)
    sched = [13, ns] if ns > 13 else [ns]
    inp = inp.model_copy(update=dict(ns_array=np.array(sched),
                                     ftol_array=np.array([1e-14, a.ftol][-len(sched):]),
                                     niter_array=np.array([6000, a.niter][-len(sched):])))
    st = settings_of(dict(M_pol=float(a.M_pol), N_tor=float(a.N_tor), lines=a.lines,
                          n_gauss=a.n_gauss, n_eta=16, n_gauss_well=16), int(inp.nfp))
    b0 = np.stack([np.asarray(inp.rbc, float), np.asarray(inp.zbs, float)])
    mask = free_mask(b0.shape, int(inp.ntor), int(a.free_m), int(a.free_m))
    vol, asp = (float(v) for v in lo.boundary_geometry(b0, int(inp.nfp), int(inp.ntor)))
    scale = 1e-3 * b0[0, 0, int(inp.ntor)] / asp
    s_list = [float(v) for v in a.s.split(",")]
    # Tight bands so that the inequality rows are active at the perturbed point.
    w = lo.Weights(iota_margin=0.0, mirror_band=0.01)
    work = tempfile.mkdtemp(prefix="check_jac_")
    p = lo.Problem(inp, s_list, st, a.mboz, a.nboz, a.objective, mask, b0, scale, w, work,
                   max_threads=a.threads)
    rng = np.random.default_rng(20261023)
    x = a.x_scale * rng.standard_normal(p.n)
    t = time.perf_counter()
    ev = p.evaluate(x, "J")
    t_J = time.perf_counter() - t
    t = time.perf_counter()
    eg = p.evaluate(x, "g")
    t_g = time.perf_counter() - t
    d = rng.standard_normal(p.n)
    d /= np.linalg.norm(d)
    Jd = ev.J @ d
    T = len(lo.OBJECTIVES[a.objective])
    n_obj = sum(len(lo.TERMS[tn][0]) if lo.TERMS[tn][1] == "lin" else 1
                for tn in lo.OBJECTIVES[a.objective]) * len(s_list)
    groups = {"objective": slice(0, n_obj), "equality": slice(n_obj, n_obj + 2),
              "iota": slice(n_obj + 2, n_obj + 2 + 2 * ns), "mirror": slice(n_obj + 2 + 2 * ns, None)}
    base = p._forward(p.boundary(x))
    f, _ = p._residual_jit("pen")
    prm = p.params()
    prm.pop("mode")

    def frozen(xx):
        e = p._forward(p.boundary(xx))
        v = base._chain.values(e._harm)
        return np.asarray(f(jnp.asarray(v), jnp.asarray(p.boundary(xx)), e._harm, prm))

    def err(fd):
        out = {}
        for g, sl in groups.items():
            ref = np.max(np.abs(Jd[sl]), initial=0.0)
            out[g] = float(np.max(np.abs(fd[sl] - Jd[sl])) / ref) if ref > 0 else None
        return out

    rows = []
    for h in (float(v) for v in a.hs.split(",")):
        fd_r = (p.evaluate(x + h * d, "r").r - p.evaluate(x - h * d, "r").r) / (2 * h)
        fd_f = (frozen(x + h * d) - frozen(x - h * d)) / (2 * h)
        rows.append(dict(h=h, refound=err(fd_r), frozen=err(fd_f)))
        print(rows[-1], flush=True)
    res = dict(input=os.path.basename(a.input), ns=ns, objective=a.objective, n=p.n, K=int(ev.r.size),
               s=s_list, t_J=t_J, t_g=t_g, info_J=ev.info, info_g=eg.info,
               grad_vs_JTr=float(np.max(np.abs(eg.g - ev.g)) / np.max(np.abs(ev.g))),
               active={g: int(np.count_nonzero(ev.r[sl])) for g, sl in groups.items()},
               vcurve=rows, boozer_route=lo.BOOZER_ROUTE)
    json.dump(res, open(a.out, "w"), indent=1, default=float)
    print(json.dumps({k: res[k] for k in ("n", "K", "t_J", "t_g", "grad_vs_JTr", "active")}))


if __name__ == "__main__":
    main()
