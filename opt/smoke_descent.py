"""A few steps of gradient descent with backtracking, as a smoke test of the driver.

The objective is the squared offset ``sum_s L_s^2``, which has a minimum, while
the plain offset ``sum_s L_s`` decreases without bound once ``Lambda_B < 0``.

The boundary is moved along the negative gradient, scaled so that the largest
coefficient change is ``step`` times the largest non-axis coefficient, and the
step is halved until the objective decreases (Armijo, c = 1e-4).  rabe's
discrete structure is recomputed at every accepted boundary.

usage: smoke_descent.py INPUT NS S1,S2,.. M_POL N_TOR MBOZ NBOZ LINES NGAUSS NSTEPS STEP OUT.json
"""

from __future__ import annotations

import json
import sys
import time

import numpy as np
from pyrabe import spectrum_derivative as sd

from check_gradient import load_input
from rabe_boundary import BoundaryObjective


def main(argv):
    path, ns, s_txt, M, N, mboz, nboz, lines, ng, nsteps, step, out = argv[:12]
    s_list = [float(v) for v in s_txt.split(",")]
    st = sd.DerivativeSettings(M_pol=float(M), N_tor=float(N), max_n_fieldlines=int(lines),
                               n_gauss=int(ng))
    obj = BoundaryObjective(load_input(path, int(ns), ftol=1e-18), s_list, st, int(mboz),
                            int(nboz), nu_star=1e-3, squared=True)
    x = obj.x0.copy()
    scale = float(np.max(np.abs(x[:, 1:, :])))
    hist = []
    J = obj.objective(x)
    for k in range(int(nsteps)):
        t = time.perf_counter()
        g = obj.grad(x)
        t_grad = time.perf_counter() - t
        p = -g / np.max(np.abs(g)) * float(step) * scale
        slope = float(np.sum(g * p))
        a, accepted = 1.0, False
        for _ in range(8):
            Jn = obj.objective(x + a * p)
            if Jn <= J + 1e-4 * a * slope:
                accepted = True
                break
            a *= 0.5
        hist.append({"iter": k, "J": J, "J_new": Jn, "alpha": a, "accepted": accepted,
                     "pred": a * slope, "t_grad": t_grad,
                     "values": obj._last.values.tolist()})
        print(f"it {k}: J {J:+.6e} -> {Jn:+.6e} (pred {a * slope:+.3e}, actual {Jn - J:+.3e}) "
              f"alpha {a:.3g} grad {t_grad:.1f}s", flush=True)
        if not accepted:
            break
        x, J = x + a * p, Jn
    with open(out, "w") as f:
        json.dump({"input": path, "ns": int(ns), "s": s_list, "history": hist,
                   "solves": obj.n_solves}, f, indent=1)


if __name__ == "__main__":
    main(sys.argv[1:])
