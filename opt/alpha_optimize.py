"""Gradient-based boundary optimization with rabe's criteria, checkpointed, with SIMPLE checks.

The experiments share the constraints and differ only in the objective
(:func:`rabe_boundary.alpha_reduce`, each term divided by its value at the start):

- ``E1`` (bootstrap): ``w1 nu*_crit + w2 (Lambda_A^2 + Lambda_B^2)``,
- ``E2`` (alpha): ``w1 nu*_crit + w2b Lambda_alpha^2 + w3 Gamma_well``,
- ``E3`` (alpha, direct): ``w1 nu*_crit + var_dJ + Gamma_s``, the label variance of
  ``J`` of the transitioning particles and the drift in the wells normalized by
  ``dJ/ds`` (needs ``ds_radial > 0``, rabe ``src/derivative/alpha_criteria.md``).

Constraints are quadratic penalties around the start point: the aspect ratio
and the volume (two-sided), ``iotaf`` inside ``[min iota0 - margin, max iota0 +
margin]`` and the Boozer mirror ratio of every objective surface inside a
relative band around its start value (:class:`rabe_boundary.MirrorRatioPenalty`).

The free variables are the boundary coefficients ``rbc``, ``zbs`` with
``m <= free_m`` and ``|n| <= free_n`` except ``rbc(0,0)`` and the unused
``m = 0, n <= 0`` entries, scaled as ``x = (b - b0) / (step_scale a_minor)``,
so that the first L-BFGS-B step (unit length in ``x``) moves the boundary by
``step_scale`` minor radii.  scipy L-BFGS-B minimizes ``F + penalties`` with the
exact gradient (one VMEC++ solve, one pyrabe VJP, one VMEC++ adjoint per
evaluation).  A failed evaluation (VMEC++ or rabe error, non-finite value)
returns a large value so that the line search backtracks.

Every evaluation is appended to ``evals.jsonl``, every accepted iteration to
``history.jsonl`` together with ``checkpoint.npz`` and the Boozer spectra of the
objective surfaces in ``spectra/iter_KKKK.npz``.  A restart continues from the
checkpoint with the start-point normalization of ``start.json`` (the L-BFGS
memory is not restored).  At the start, every ``simple_every``-th iteration and
at the end, the boundary is re-solved, its wout written to
``simple/<tag>/wout.nc`` and one SIMPLE alpha-loss run per ``simple_s`` staged
and submitted with ``sbatch <simple_sbatch> <rundir>`` (if configured).

usage:
  alpha_optimize.py run CONFIG.json
  alpha_optimize.py ftol-check CONFIG.json FTOL1,FTOL2,..   (writes ftol_choice.json)
  alpha_optimize.py finalize CONFIG.json                    (SIMPLE for the last checkpoint)
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
import traceback

import jax.numpy as jnp
import numpy as np
import scipy.optimize
import vmecpp
from pyrabe import spectrum_derivative as sd

from rabe_boundary import (AlphaWeights, BoundaryObjective, MirrorRatioPenalty, alpha_reduce,
                           alpha_terms, aspect_penalty, boundary_geometry, iota_range_penalty,
                           volume_penalty)

EXPERIMENTS = {
    "E1": dict(nu_star_crit=1.0, offset=1.0, Lambda_alpha=0.0, Gamma_well=0.0),
    "E2": dict(nu_star_crit=1.0, offset=0.0, Lambda_alpha=1.0, Gamma_well=1.0),
    "E3": dict(nu_star_crit=1.0, offset=0.0, Lambda_alpha=0.0, Gamma_well=0.0, var_dJ=1.0,
               Gamma_s=1.0),
}

DEFAULTS = dict(
    ns=75, ftol=1e-18, niter=200000, ns_coarse=13, niter_coarse=6000,
    s_list=[0.25, 0.4, 0.55, 0.7], M_pol=0.0, N_tor=None, mboz=17, nboz=16,
    lines=100, n_gauss=64, n_eta=32, n_gauss_well=32, ds_radial=0.0,
    free_m=4, free_n=4, step_scale=1e-3, bound=50.0, maxiter=60, threads=8,
    w_aspect=1e3, w_volume=1e3, iota_margin=0.02, w_iota=1e3,
    mirror_band=0.1, w_mirror=1e3,
    simple_every=5, simple_s=[0.25, 0.5], simple_n=5000, simple_time=0.2, simple_sbatch=None,
)

SIMPLE_TEMPLATE = """&config
notrace_passing = 0
nper = 1000
npoiper = 100
ntimstep = 10000
ntestpart = {n}
trace_time = {tt}d0
sbeg = {sbeg}d0
num_surf = 1
phibeg = 0.d0
thetabeg = 0.d0
contr_pp = -1.0d0
facE_al = 1.0
npoiper2 = 128
n_e = 2
n_d = 4
netcdffile = 'wout.nc'
ns_s = 5
ns_tp = 5
multharm = 5
isw_field_type = 2
generate_start_only = .False.
startmode = 1
grid_density = 0d0
integmode = 1
relerr = 1d-13
tcut = -1d0
fast_class = .False.
vmec_B_scale = {b:.15e}
vmec_RZ_scale = {rz:.15e}
swcoll = .False.
deterministic = .True.
old_axis_healing = .True.
old_axis_healing_boundary = .True.
batch_size = 30000
ran_seed = 12345
/
"""
A_TARGET, B_TARGET = 1.704, 5.865  # ARIES-CS minor radius [m] and <B> [T]


def load_config(path):
    cfg = dict(DEFAULTS)
    cfg.update(json.load(open(path)))
    cfg["outdir"] = os.path.abspath(cfg.get("outdir", os.path.dirname(os.path.abspath(path))))
    if cfg["ftol"] == "auto":
        choice = os.path.join(cfg["outdir"], cfg.get("ftol_choice", "../ftol_check/ftol_choice.json"))
        cfg["ftol"] = float(json.load(open(choice))["ftol"])
    return cfg


def load_input(cfg, ftol=None):
    inp = vmecpp.VmecInput.from_file(cfg["input"])
    ftol = cfg["ftol"] if ftol is None else ftol
    ns = [cfg["ns_coarse"], cfg["ns"]] if cfg["ns_coarse"] else [cfg["ns"]]
    ft = [1e-14, ftol] if cfg["ns_coarse"] else [ftol]
    it = [cfg["niter_coarse"], cfg["niter"]] if cfg["ns_coarse"] else [cfg["niter"]]
    return inp.model_copy(update={"ns_array": np.array(ns), "ftol_array": np.array(ft),
                                  "niter_array": np.array(it),
                                  "return_outputs_even_if_not_converged": True})


def settings_of(cfg, nfp):
    N = float(nfp if cfg["N_tor"] is None else cfg["N_tor"])
    return sd.DerivativeSettings(M_pol=float(cfg["M_pol"]), N_tor=N,
                                 max_n_fieldlines=int(cfg["lines"]), n_gauss=int(cfg["n_gauss"]),
                                 n_eta=int(cfg["n_eta"]), n_gauss_well=int(cfg["n_gauss_well"]))


def free_mask(shape, ntor, free_m, free_n):
    _, mpol, nn = shape
    m = np.arange(mpol)[:, None] * np.ones((1, nn), int)
    n = (np.arange(nn) - ntor)[None, :] * np.ones((mpol, 1), int)
    one = (m <= free_m) & (np.abs(n) <= free_n) & ~((m == 0) & (n <= 0))
    mask = np.stack([one, one.copy()])
    mask[0, 0, ntor] = False  # rbc(0,0) stays: the size is held by aspect and volume
    return mask


def json_line(path, rec):
    with open(path, "a") as f:
        f.write(json.dumps(rec) + "\n")


class Run:
    def __init__(self, cfg, ftol=None):
        self.cfg = cfg
        os.makedirs(cfg["outdir"], exist_ok=True)
        self.inp = load_input(cfg, ftol)
        self.nfp, self.ntor = int(self.inp.nfp), int(self.inp.ntor)
        self.st = settings_of(cfg, self.nfp)
        work = os.path.join(cfg["outdir"], "work")
        os.makedirs(work, exist_ok=True)
        self.obj = BoundaryObjective(self.inp, cfg["s_list"], self.st, cfg["mboz"], cfg["nboz"],
                                     workdir=work, max_threads=int(cfg["threads"]),
                                     ds_radial=float(cfg["ds_radial"]))
        self.b0 = self.obj.x0
        vol, asp = (float(v) for v in boundary_geometry(self.b0, self.nfp, self.ntor))
        rmaj = float(self.b0[0, 0, self.ntor])
        self.a_minor = rmaj / asp
        self.mask = free_mask(self.b0.shape, self.ntor, cfg["free_m"], cfg["free_n"])
        self.scale = cfg["step_scale"] * self.a_minor
        self.mirror = MirrorRatioPenalty(cfg["s_list"], 0.0, np.inf, 0.0, cfg["mboz"],
                                         cfg["nboz"]).bind(self.obj.base_wout)
        self.vol0, self.asp0 = vol, asp

    # -- start point, normalization and constraint targets -------------------------------
    def start(self, experiment):
        path = os.path.join(self.cfg["outdir"], "start.json")
        if os.path.exists(path):
            st = json.load(open(path))
        else:
            ev = self.obj.evaluate(self.b0)
            iota0 = np.asarray(ev.harm["iotaf"])
            d0 = np.asarray(self.mirror.ratios(ev.harm))
            terms0 = np.asarray(alpha_terms(jnp.asarray(ev.values)))
            st = dict(values=ev.values.tolist(), terms=terms0.tolist(), aspect=self.asp0,
                      volume=self.vol0, iota_min=float(iota0.min()), iota_max=float(iota0.max()),
                      mirror=d0.tolist(), a_minor=self.a_minor, solver=self.obj.solver_info)
            json.dump(st, open(path, "w"), indent=1)
        c = self.cfg
        if EXPERIMENTS[experiment].get("Gamma_s", 0.0) > 0.0 and not c["ds_radial"] > 0.0:
            raise ValueError(f"{experiment} uses Gamma_s and needs ds_radial > 0")
        if np.asarray(st["values"]).shape[1] != len(sd.OUTPUT_NAMES):
            raise ValueError("start.json was written by an older rabe/pyrabe (other outputs)")
        w = AlphaWeights(**EXPERIMENTS[experiment]).normalized(np.asarray(st["values"]))
        d0 = np.asarray(st["mirror"])
        self.mirror.lower, self.mirror.upper = d0 * (1 - c["mirror_band"]), d0 * (1 + c["mirror_band"])
        self.mirror.weight = c["w_mirror"] / float(np.mean(d0)) ** 2
        self.named = {
            "aspect": aspect_penalty(st["aspect"], c["w_aspect"] / st["aspect"] ** 2, self.nfp,
                                     self.ntor, upper_only=False),
            "volume": volume_penalty(st["volume"], c["w_volume"], self.nfp, self.ntor),
            "iota": iota_range_penalty(st["iota_min"] - c["iota_margin"],
                                       st["iota_max"] + c["iota_margin"], c["w_iota"]),
            "mirror": self.mirror,
        }
        self.obj.set_objective(alpha_reduce(w), list(self.named.values()))
        self.weights = w
        return st

    def boundary(self, x):
        b = self.b0.copy()
        b[self.mask] += self.scale * np.asarray(x)
        return b

    def record(self, b, F, ok, err=None):
        ev = self.obj._last
        rec = dict(t=time.time(), ok=ok, F=F, n_solves=self.obj.n_solves, err=err)
        if ok and ev is not None and np.array_equal(ev.boundary, b):
            terms = np.asarray(alpha_terms(jnp.asarray(ev.values)))
            pen = {k: float(c(jnp.asarray(b), ev.harm)) for k, c in self.named.items()}
            vol, asp = (float(v) for v in boundary_geometry(b, self.nfp, self.ntor))
            iota = np.asarray(ev.harm["iotaf"])
            rec.update(values=ev.values.tolist(), terms=terms.tolist(),
                       objective=float(self.obj.reduce(jnp.asarray(ev.values))), penalties=pen,
                       aspect=asp, volume=vol, iota_min=float(iota.min()),
                       iota_max=float(iota.max()),
                       mirror=np.asarray(self.mirror.ratios(ev.harm)).tolist(),
                       solver=self.obj.solver_info, timing=ev.timing)
        return rec

    # -- SIMPLE ------------------------------------------------------------------------
    def simple(self, b, tag):
        c = self.cfg
        root = os.path.join(c["outdir"], "simple", tag)
        os.makedirs(root, exist_ok=True)
        wout_path = os.path.join(root, "wout.nc")
        if not os.path.exists(wout_path):
            out = self.obj.solve(b)
            out.wout.save(wout_path)
        import netCDF4
        with netCDF4.Dataset(wout_path) as d:
            rz = A_TARGET / float(np.asarray(d["Aminor_p"][...]))
            bs = B_TARGET / abs(float(np.asarray(d["volavgB"][...])))
            conv = {k: float(np.asarray(d[k][...])) for k in ("fsqr", "fsqz", "fsql")}
        jobs = []
        for s in c["simple_s"]:
            run = os.path.join(root, f"s{s}_n{c['simple_n']}")
            os.makedirs(run, exist_ok=True)
            if not os.path.exists(os.path.join(run, "wout.nc")):
                os.symlink("../wout.nc", os.path.join(run, "wout.nc"))
            open(os.path.join(run, "simple.in"), "w").write(SIMPLE_TEMPLATE.format(
                n=c["simple_n"], tt=c["simple_time"], sbeg=s, b=bs, rz=rz))
            job = None
            if c["simple_sbatch"] and not os.path.exists(os.path.join(run, "submitted")):
                r = subprocess.run(["sbatch", "--parsable", c["simple_sbatch"], run],
                                   capture_output=True, text=True)
                job = r.stdout.strip() or r.stderr.strip()
                open(os.path.join(run, "submitted"), "w").write(job + "\n")
            jobs.append(dict(run=run, job=job))
        json_line(os.path.join(c["outdir"], "simple_jobs.jsonl"),
                  dict(tag=tag, t=time.time(), wout_conv=conv, jobs=jobs))
        return jobs


def cmd_run(cfg):
    run = Run(cfg)
    st = run.start(cfg["experiment"])
    out = cfg["outdir"]
    os.makedirs(os.path.join(out, "spectra"), exist_ok=True)
    ck = os.path.join(out, "checkpoint.npz")
    x0 = np.zeros(int(run.mask.sum()))
    k0 = 0
    if os.path.exists(ck):
        z = np.load(ck)
        x0, k0 = np.asarray(z["x"]), int(z["iteration"])
        print(f"restart from iteration {k0}", flush=True)
    else:
        run.simple(run.b0, "iter0000")
    # F = number of active normalized terms at the start; a failed evaluation returns 10x that.
    F_fail = 10.0 * max(1.0, sum(v > 0 for v in EXPERIMENTS[cfg["experiment"]].values()))
    state = dict(k=k0, g_last=np.zeros_like(x0), last=None)

    def fun(x):
        b = run.boundary(x)
        t = time.perf_counter()
        try:
            F = run.obj.objective(b)
            g = run.obj.grad(b)[run.mask] * run.scale
            if not (np.isfinite(F) and np.all(np.isfinite(g))):
                raise FloatingPointError("non-finite objective or gradient")
            rec = run.record(b, F, True)
            state["g_last"], state["last"] = g, (x.copy(), b, F, rec)
        except Exception as e:  # noqa: BLE001 - a failed evaluation makes the line search backtrack
            rec = run.record(b, F_fail, False, f"{type(e).__name__}: {e}")
            traceback.print_exc()
            F, g = F_fail, state["g_last"]
        rec["wall"] = time.perf_counter() - t
        rec["iteration"] = state["k"]
        json_line(os.path.join(out, "evals.jsonl"), rec)
        print(f"eval it {state['k']}: F {F:.6e} ok {rec['ok']} wall {rec['wall']:.0f}s", flush=True)
        return F, g

    def callback(intermediate_result):
        x = np.asarray(intermediate_result.x)
        k = state["k"]
        last = state["last"]
        if last is None or not np.array_equal(last[0], x):
            fun(x)  # make the stored evaluation that of the accepted point
            last = state["last"]
        _, b, F, rec = last
        np.savez(ck, x=x, boundary=b, iteration=k, F=F, mask=run.mask, b0=run.b0,
                 scale=run.scale)
        ev = run.obj._last
        spec = {}
        for i in range(len(cfg["s_list"])):
            p, xm, xn = run.obj._last.chain._params(i, ev.harm)
            spec[f"bmnc_{i}"], spec["xm"], spec["xn"] = (np.asarray(p.bmnc), np.asarray(xm),
                                                         np.asarray(xn))
        np.savez(os.path.join(out, "spectra", f"iter_{k:04d}.npz"), s=np.asarray(cfg["s_list"]),
                 boundary=b, **spec)
        json_line(os.path.join(out, "history.jsonl"), dict(rec, iteration=k))
        print(f"iteration {k}: F {F:.6e}", flush=True)
        if cfg["simple_every"] and k > 0 and k % int(cfg["simple_every"]) == 0:
            run.simple(b, f"iter{k:04d}")
        state["k"] = k + 1

    # Iteration 0 record (start point or restart point) before the first step.
    callback(scipy.optimize.OptimizeResult(x=x0, fun=None))
    res = scipy.optimize.minimize(fun, x0, jac=True, method="L-BFGS-B", callback=callback,
                                  bounds=[(-cfg["bound"], cfg["bound"])] * x0.size,
                                  options=dict(maxiter=max(int(cfg["maxiter"]) - k0, 1), maxcor=20,
                                               ftol=1e-10, gtol=1e-8))
    json.dump(dict(message=str(res.message), success=bool(res.success), nit=int(res.nit),
                   nfev=int(res.nfev), F=float(res.fun), n_solves=run.obj.n_solves),
              open(os.path.join(out, "result.json"), "w"), indent=1)
    run.simple(run.boundary(res.x), "final")


def cmd_finalize(cfg):
    """SIMPLE for the last checkpoint, for example after a time limit or a crash."""
    if os.path.exists(os.path.join(cfg["outdir"], "simple", "final")):
        print("final SIMPLE already staged")
        return
    z = np.load(os.path.join(cfg["outdir"], "checkpoint.npz"))
    run = Run(cfg)
    run.simple(np.asarray(z["boundary"]), f"last{int(z['iteration']):04d}")


def cmd_ftol_check(cfg, ftols):
    """Values and E1/E2 gradients at the start boundary for several ftol.

    The reference is the smallest ftol.  The chosen ftol is the largest one whose
    E1 and E2 gradients agree with the reference to ``1e-2`` (relative 2-norm) and
    whose terms agree to ``1e-3`` (relative, per surface and term).
    """
    out = cfg["outdir"]
    os.makedirs(out, exist_ok=True)
    rows = []
    for ftol in sorted(ftols):
        t = time.perf_counter()
        run = Run(dict(cfg, outdir=os.path.join(out, f"ftol_{ftol:.0e}")), ftol=ftol)
        row = dict(ftol=ftol)
        for e in ("E1", "E2"):
            st = run.start(e)
            g = run.obj.grad(run.b0)[run.mask] * run.scale
            row[e] = dict(F=run.obj.objective(run.b0), g=g.tolist())
        row.update(terms=st["terms"], values=st["values"], solver=st["solver"],
                   wall=time.perf_counter() - t)
        rows.append(row)
        json.dump(rows, open(os.path.join(out, "ftol_check.json"), "w"), indent=1)
        print(f"ftol {ftol:.0e}: solver {st['solver']} wall {row['wall']:.0f}s", flush=True)
    ref = rows[0]
    choice = ref["ftol"]
    for row in rows:
        dg = max(np.linalg.norm(np.subtract(row[e]["g"], ref[e]["g"]))
                 / np.linalg.norm(ref[e]["g"]) for e in ("E1", "E2"))
        dt = float(np.max(np.abs(np.subtract(row["terms"], ref["terms"]))
                          / np.maximum(np.abs(ref["terms"]), 1e-300)))
        row["dgrad"], row["dterms"] = float(dg), dt
        if dg <= 1e-2 and dt <= 1e-3:
            choice = max(choice, row["ftol"])
    choice = min(choice, 1e-16)
    json.dump(rows, open(os.path.join(out, "ftol_check.json"), "w"), indent=1)
    json.dump(dict(ftol=choice, rule="largest ftol with |dg|/|g| <= 1e-2 (E1, E2) and terms "
                   "<= 1e-3 relative vs the smallest ftol, capped at 1e-16",
                   table=[{k: r[k] for k in ("ftol", "dgrad", "dterms", "wall")} for r in rows]),
              open(os.path.join(out, "ftol_choice.json"), "w"), indent=1)
    print("ftol choice", choice, flush=True)


def main(argv):
    cmd, path = argv[0], argv[1]
    cfg = load_config(path) if cmd != "ftol-check" else dict(DEFAULTS, **json.load(open(path)))
    if cmd == "ftol-check":
        cfg["outdir"] = os.path.abspath(cfg.get("outdir", os.path.dirname(os.path.abspath(path))))
        cmd_ftol_check(cfg, [float(v) for v in argv[2].split(",")])
    elif cmd == "run":
        cmd_run(cfg)
    elif cmd == "finalize":
        cmd_finalize(cfg)
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
