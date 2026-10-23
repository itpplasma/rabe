"""Campaign entry point: derivative-based boundary optimization with rabe's criteria.

usage:
  run_opt.py --engine vmecpp --objective E1|E2|E3 --solver lm|trf|lbfgs|al \\
             --input INPUT --outdir DIR [--config CONFIG.json] [options]

Stages continue in the free boundary modes and in ``ns`` (``--stages
"m:n:ns:maxiter,..."``, default ``2:2:25:15,3:3:51:15,4:4:75:20,6:6:75:20``):
each stage frees ``rbc, zbs`` with ``m <= m_s``, ``|n| <= n_s`` and solves VMEC++
at ``ns_s`` (cold from ``13 -> ns_s`` the first time, then hot-restarted from
the previous solve), starting from the boundary the previous stage ended at.
Each ``ns`` has its own start normalization and constraint targets, taken at
the start boundary (``reference_ns<N>.json``).  Every evaluation goes to
``evals.jsonl``, every accepted iteration to ``history.jsonl`` and
``checkpoint.json`` (boundary, solver state); a rerun with the same
``--outdir`` resumes at the checkpoint.  ``--start-index k --perturb p`` starts
from a random perturbation of norm ``~p sqrt(n)`` in ``x`` (multi-start, one
directory per start).  ``--simple-sbatch`` stages and submits SIMPLE for the
final boundary (and the start with ``--simple-start true``) like ``alpha_optimize.py``.

Boozer route of all criteria: pyrabe.boozer_libneo (libneo hybrid transform of
the VMEC++ primaries); recorded in ``meta.json``.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np
import scipy.optimize
import vmecpp

import lsq_opt as lo
import vmecpp_jacobian as vj
from alpha_optimize import A_TARGET, B_TARGET, SIMPLE_TEMPLATE, free_mask, settings_of

DEFAULTS = dict(
    engine="vmecpp", objective="E1", solver="lm",
    stages="2:2:25:15,3:3:51:15,4:4:75:20,6:6:75:20",
    ftol=1e-18, niter=200000, ns_coarse=13, niter_coarse=6000,
    s_list=[0.25, 0.4, 0.55, 0.7], M_pol=0.0, N_tor=None, mboz=17, nboz=16,
    lines=100, n_gauss=64, n_eta=32, n_gauss_well=32,
    step_scale=1e-3, bound=50.0, max_step=5.0, lam0=1e-2, threads=8,
    w_aspect=1e3, w_volume=1e3, w_iota=1e3, w_mirror=1e3, iota_margin=0.02, mirror_band=0.1,
    al_outer=6, al_inner=6, al_mu0=1e3, al_tol=1e-3,
    start_index=0, perturb=0.0, seed=20261023, warm=True,
    simple_sbatch=None, simple_start=False, simple_s=[0.25, 0.5], simple_n=5000, simple_time=0.2,
)


def parse(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config")
    for k, v in DEFAULTS.items():
        kw = dict(default=None)
        if isinstance(v, bool):
            kw["type"] = lambda s: s.lower() in ("1", "true", "yes")
        elif isinstance(v, list):
            kw["type"] = lambda s: [float(t) for t in s.split(",")]
        elif isinstance(v, (int, float)) and not isinstance(v, bool):
            kw["type"] = float if isinstance(v, float) else int
        ap.add_argument("--" + k.replace("_", "-"), dest=k, **kw)
    ap.add_argument("--input", required=False)
    ap.add_argument("--outdir", required=False)
    a = ap.parse_args(argv)
    cfg = dict(DEFAULTS)
    if a.config:
        cfg.update(json.load(open(a.config)))
    cfg.update({k: v for k, v in vars(a).items() if v is not None and k != "config"})
    if not cfg.get("input") or not cfg.get("outdir"):
        ap.error("--input and --outdir (or both in --config) are required")
    if cfg["engine"] != "vmecpp":
        raise SystemExit(f"engine {cfg['engine']} is not wired; only vmecpp")
    if cfg["objective"] not in lo.OBJECTIVES or cfg["solver"] not in ("lm", "trf", "lbfgs", "al"):
        ap.error("objective E1|E2|E3, solver lm|trf|lbfgs|al")
    cfg["outdir"] = os.path.abspath(cfg["outdir"])
    if cfg["ftol"] == "auto" or cfg.get("ftol_choice"):
        path = cfg.get("ftol_choice")
        cfg["ftol"] = float(json.load(open(path))["ftol"])
    return cfg


def stages_of(spec):
    out = []
    for item in str(spec).split(","):
        m, n, ns, it = (int(v) for v in item.split(":"))
        out.append(dict(m=m, n=n, ns=ns, maxiter=it))
    return out


def padded_input(inp, mpol, ntor):
    """``inp`` with at least ``mpol`` poloidal and ``ntor`` toroidal modes (zero padded)."""
    if inp.mpol >= mpol and inp.ntor >= ntor:
        return inp
    M, N = max(inp.mpol, mpol), max(inp.ntor, ntor)
    upd = {"mpol": M, "ntor": N}
    for name in ("rbc", "zbs", "rbs", "zbc"):
        a = getattr(inp, name)
        if a is None:
            continue
        a = np.asarray(a)
        b = np.zeros((M, 2 * N + 1))
        b[:a.shape[0], N - inp.ntor:N + inp.ntor + 1] = a
        upd[name] = b
    for name in ("raxis_c", "zaxis_s", "raxis_s", "zaxis_c"):
        a = getattr(inp, name)
        if a is not None:
            b = np.zeros(N + 1)
            b[:np.size(a)] = np.asarray(a)
            upd[name] = b
    return inp.model_copy(update=upd)


def stage_input(base, cfg, ns):
    if cfg["ns_coarse"] and ns > cfg["ns_coarse"]:
        ns_a, ft, it = [cfg["ns_coarse"], ns], [1e-14, cfg["ftol"]], [cfg["niter_coarse"], cfg["niter"]]
    else:
        ns_a, ft, it = [ns], [cfg["ftol"]], [cfg["niter"]]
    return base.model_copy(update={"ns_array": np.array(ns_a), "ftol_array": np.array(ft),
                                   "niter_array": np.array(it)})


def jline(path, rec):
    with open(path, "a") as f:
        f.write(json.dumps(rec, default=float) + "\n")


def stage_simple(cfg, base, boundary, tag):
    """Solve ``boundary``, write its wout and SIMPLE runs; submit if configured."""
    import netCDF4
    root = os.path.join(cfg["outdir"], "simple", tag)
    os.makedirs(root, exist_ok=True)
    wout = os.path.join(root, "wout.nc")
    ns = max(s["ns"] for s in stages_of(cfg["stages"]))
    if not os.path.exists(wout):
        inp = stage_input(base, cfg, ns).model_copy(update={"rbc": boundary[0], "zbs": boundary[1]})
        vmecpp.run(inp, max_threads=int(cfg["threads"]), verbose=False).wout.save(wout)
    with netCDF4.Dataset(wout) as d:
        rz = A_TARGET / float(np.asarray(d["Aminor_p"][...]))
        bs = B_TARGET / abs(float(np.asarray(d["volavgB"][...])))
    jobs = []
    for s in cfg["simple_s"]:
        run = os.path.join(root, f"s{s}_n{cfg['simple_n']}")
        os.makedirs(run, exist_ok=True)
        if not os.path.exists(os.path.join(run, "wout.nc")):
            os.symlink("../wout.nc", os.path.join(run, "wout.nc"))
        open(os.path.join(run, "simple.in"), "w").write(SIMPLE_TEMPLATE.format(
            n=cfg["simple_n"], tt=cfg["simple_time"], sbeg=s, b=bs, rz=rz))
        job = None
        if cfg["simple_sbatch"] and not os.path.exists(os.path.join(run, "submitted")):
            r = subprocess.run(["sbatch", "--parsable", cfg["simple_sbatch"], run],
                               capture_output=True, text=True)
            job = r.stdout.strip() or r.stderr.strip()
            open(os.path.join(run, "submitted"), "w").write(job + "\n")
        jobs.append(dict(run=run, job=job))
    jline(os.path.join(cfg["outdir"], "simple_jobs.jsonl"), dict(tag=tag, jobs=jobs, t=time.time()))


class Campaign:
    def __init__(self, cfg):
        self.cfg = cfg
        self.out = cfg["outdir"]
        os.makedirs(self.out, exist_ok=True)
        self.stages = stages_of(cfg["stages"])
        inp = vmecpp.VmecInput.from_file(cfg["input"])
        mmax = max(s["m"] for s in self.stages) + 1
        nmax = max(s["n"] for s in self.stages)
        self.base = padded_input(inp, mmax, nmax)
        self.nfp, self.ntor = int(self.base.nfp), int(self.base.ntor)
        self.b0 = np.stack([np.asarray(self.base.rbc, float), np.asarray(self.base.zbs, float)])
        vol, asp = (float(v) for v in lo.boundary_geometry(self.b0, self.nfp, self.ntor))
        self.scale = cfg["step_scale"] * float(self.b0[0, 0, self.ntor]) / asp
        self.settings = settings_of(cfg, self.nfp)
        self.weights = lo.Weights(cfg["w_aspect"], cfg["w_volume"], cfg["w_iota"], cfg["w_mirror"],
                                  cfg["iota_margin"], cfg["mirror_band"])
        self.ck_path = os.path.join(self.out, "checkpoint.json")
        self.ck = json.load(open(self.ck_path)) if os.path.exists(self.ck_path) else None
        json.dump(dict(cfg=cfg, boozer_route=lo.BOOZER_ROUTE, vmecpp=vj._version,
                       mpol=int(self.base.mpol), ntor=self.ntor, scale=self.scale,
                       output_names=list(lo.sd.OUTPUT_NAMES)),
                  open(os.path.join(self.out, "meta.json"), "w"), indent=1, default=float)

    def problem(self, k):
        st, cfg = self.stages[k], self.cfg
        ref_path = os.path.join(self.out, f"reference_ns{st['ns']}.json")
        ref = lo.Reference.from_json(json.load(open(ref_path))) if os.path.exists(ref_path) else None
        mask = free_mask(self.b0.shape, self.ntor, st["m"], st["n"])
        p = lo.Problem(stage_input(self.base, cfg, st["ns"]), cfg["s_list"], self.settings,
                       cfg["mboz"], cfg["nboz"], cfg["objective"], mask, self.b0, self.scale,
                       self.weights, os.path.join(self.out, "work"), reference=ref,
                       warm_key=(self.out, st["ns"]) if cfg["warm"] else None,
                       max_threads=int(cfg["threads"]))
        if ref is None:
            json.dump(p.reference.to_json(), open(ref_path, "w"), indent=1)
        return p

    def save(self, k, boundary, solver_state, ev):
        ck = dict(stage=k, boundary=np.asarray(boundary).tolist(), solver_state=solver_state,
                  F=ev.F, objective=ev.objective, violation=ev.violation(), t=time.time())
        tmp = self.ck_path + ".tmp"
        json.dump(ck, open(tmp, "w"), default=lambda o: np.asarray(o).tolist())
        os.replace(tmp, self.ck_path)
        self.ck = ck

    def run(self):
        cfg = self.cfg
        k0, boundary, sstate = 0, self.b0.copy(), None
        if self.ck:
            k0, boundary = int(self.ck["stage"]), np.asarray(self.ck["boundary"])
            sstate = self.ck.get("solver_state")
            if sstate and sstate.get("stage_done"):
                k0, sstate = k0 + 1, None
            print(f"resume at stage {k0}", flush=True)
        elif cfg["simple_sbatch"] and cfg["simple_start"]:
            stage_simple(cfg, self.base, self.b0, "start")
        for k in range(k0, len(self.stages)):
            boundary = self.run_stage(k, boundary, sstate)
            sstate = None
        json.dump(dict(final_boundary=boundary.tolist(), stats=vj.STATS.as_dict()),
                  open(os.path.join(self.out, "result.json"), "w"), indent=1, default=float)
        if cfg["simple_sbatch"]:
            stage_simple(cfg, self.base, boundary, "final")
        return boundary

    def run_stage(self, k, boundary, sstate):
        cfg, st = self.cfg, self.stages[k]
        t_stage = time.perf_counter()
        p = self.problem(k)
        x0 = p.x_of(boundary)
        if k == 0 and sstate is None and not self.ck and cfg["perturb"] > 0:
            rng = np.random.default_rng([int(cfg["seed"]), int(cfg["start_index"])])
            x0 = x0 + cfg["perturb"] * rng.standard_normal(x0.size)
        evals, hist = os.path.join(self.out, "evals.jsonl"), os.path.join(self.out, "history.jsonl")
        counters = dict(nfev=0, t0=time.perf_counter())

        def on_eval(ev, kind):
            counters["nfev"] += 1
            jline(evals, dict(stage=k, kind=kind, n=p.n, wall=time.perf_counter() - counters["t0"],
                              **lo.summary(ev)))

        def on_iter(ev, s):
            state = {kk: v for kk, v in s.items() if kk in
                     ("lam", "nu", "k", "outer", "mu", "lam_eq", "lam_in", "viol", "lm")}
            state["x"] = np.asarray(ev.x).tolist()
            self.save(k, ev.boundary, state, ev)
            jline(hist, dict(stage=k, iteration=s.get("k", s.get("outer")), al=bool(s.get("al")),
                             nfev=counters["nfev"], wall=time.perf_counter() - counters["t0"],
                             F=ev.F, objective=ev.objective, violation=ev.violation(),
                             terms=ev.terms.tolist(), c_eq=ev.c_eq.tolist(), lam=s.get("lam"),
                             xnorm=float(np.linalg.norm(ev.x)), stats=vj.STATS.as_dict()))
            print(f"stage {k} it {s.get('k', s.get('outer'))}: F {ev.F:.6e} obj {ev.objective:.4e} "
                  f"viol {ev.violation():.2e} nfev {counters['nfev']}", flush=True)

        lm_kw = dict(max_step=cfg["max_step"], bound=cfg["bound"], lam0=cfg["lam0"])
        restore = None
        if sstate and "x" in sstate:
            x0 = np.asarray(sstate["x"])
            restore = {kk: v for kk, v in sstate.items() if kk != "x"}
        solver = cfg["solver"]
        if solver == "lm":
            lm_state = {kk: restore[kk] for kk in ("lam", "nu", "k") if restore and kk in restore}
            ev, s, msg = lo.levenberg_marquardt(p, x0, maxiter=st["maxiter"], on_iter=on_iter,
                                                on_eval=on_eval, state=lm_state or None, **lm_kw)
        elif solver == "al":
            al_state = None
            if restore and "outer" in restore:
                al_state = {kk: restore[kk] for kk in ("outer", "mu", "lam_eq", "lam_in", "viol", "lm")
                            if kk in restore}
                al_state["x"] = x0
            ev, s, msg = lo.augmented_lagrangian(p, x0, outer=int(cfg["al_outer"]),
                                                 inner=int(cfg["al_inner"]), mu0=cfg["al_mu0"],
                                                 tol=cfg["al_tol"], on_iter=on_iter, on_eval=on_eval,
                                                 state=al_state, **lm_kw)
        else:
            ev, msg = self.scipy_stage(p, x0, st["maxiter"], solver, on_iter, on_eval)
        self.save(k, ev.boundary, dict(stage_done=True, message=msg, x=np.asarray(ev.x).tolist()), ev)
        jline(os.path.join(self.out, "stages.jsonl"),
              dict(stage=k, **st, n_free=p.n, message=msg, F=ev.F, objective=ev.objective,
                   violation=ev.violation(), wall=time.perf_counter() - t_stage,
                   n_eval=p.n_eval, stats=vj.STATS.as_dict()))
        print(f"stage {k} done: {msg}, F {ev.F:.6e}", flush=True)
        return ev.boundary

    def scipy_stage(self, p, x0, maxiter, solver, on_iter, on_eval):
        cfg = self.cfg
        best = dict(ev=None)
        it = dict(k=0)

        def keep(ev):
            if ev.ok and (best["ev"] is None or ev.F < best["ev"].F):
                best["ev"] = ev
                on_iter(ev, dict(k=it["k"]))
                it["k"] += 1

        ev0 = p.evaluate(x0, "J")
        on_eval(ev0, "J")
        if not ev0.ok:
            raise RuntimeError(f"start failed: {ev0.info.get('error')}")
        keep(ev0)
        r_fail = 10.0 * np.sqrt(max(ev0.F, 1.0) / ev0.r.size) * np.ones_like(ev0.r)
        cache = dict(J=ev0.J, g=ev0.g)
        bounds = (-cfg["bound"], cfg["bound"])
        if solver == "trf":
            def fun(x):
                ev = p.evaluate(x, "r")
                on_eval(ev, "r")
                keep(ev)
                return ev.r if ev.ok else r_fail

            def jac(x):
                ev = p.evaluate(x, "J")
                on_eval(ev, "J")
                if ev.ok:
                    cache["J"] = ev.J
                return cache["J"]
            res = scipy.optimize.least_squares(fun, x0, jac=jac, method="trf", bounds=bounds,
                                               x_scale=1.0, max_nfev=int(maxiter) * 3,
                                               tr_solver="exact", ftol=1e-10, xtol=1e-10, gtol=1e-12)
        else:
            def fg(x):
                ev = p.evaluate(x, "g")
                on_eval(ev, "g")
                keep(ev)
                if not ev.ok:
                    return 10.0 * max(ev0.F, 1.0), cache["g"]
                cache["g"] = ev.g
                return ev.F, ev.g
            res = scipy.optimize.minimize(fg, x0, jac=True, method="L-BFGS-B",
                                          bounds=[bounds] * x0.size,
                                          options=dict(maxiter=int(maxiter), maxcor=20, ftol=1e-12,
                                                       gtol=1e-10))
        return best["ev"], str(res.message)


def main(argv):
    cfg = parse(argv)
    Campaign(cfg).run()


if __name__ == "__main__":
    main(sys.argv[1:])
