"""Radial resolution scan: values, boundary gradients and cost per gradient.

usage: ns_scan.py INPUT NS1,NS2,.. S1,S2,.. M_POL N_TOR MBOZ NBOZ LINES NGAUSS OUT.json [FTOL]
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
    path, ns_txt, s_txt, M, N, mboz, nboz, lines, ng, out = argv[:10]
    ftol = float(argv[10]) if len(argv) > 10 else 1e-20
    s_list = [float(v) for v in s_txt.split(",")]
    st = sd.DerivativeSettings(M_pol=float(M), N_tor=float(N), max_n_fieldlines=int(lines),
                               n_gauss=int(ng))
    rows = []
    for ns in (int(v) for v in ns_txt.split(",")):
        obj = BoundaryObjective(load_input(path, ns, ftol=ftol), s_list, st, int(mboz), int(nboz),
                                nu_star=1e-3)
        x0 = obj.x0
        t = time.perf_counter()
        J = obj.objective(x0)
        g = obj.grad(x0)
        dt = time.perf_counter() - t
        ev = obj._last
        rows.append({"ns": ns, "J": J, "values": ev.values.tolist(), "grad": g.tolist(),
                     "timing": dict(ev.timing, value_and_grad=dt)})
        print(f"ns={ns:4d} J={J:+.6e} t={dt:6.1f}s "
              + " ".join(f"{k}={v:.1f}" for k, v in ev.timing.items())
              + "  LA=" + ",".join(f"{v:+.4e}" for v in ev.values[:, 0])
              + "  LB=" + ",".join(f"{v:+.4e}" for v in ev.values[:, 1]), flush=True)
        with open(out, "w") as f:
            json.dump({"input": path, "s": s_list, "ftol": ftol, "rows": rows}, f)
    ref = np.asarray(rows[-1]["grad"])
    for r in rows:
        g = np.asarray(r["grad"])
        cos = float(np.sum(g * ref) / np.linalg.norm(g) / np.linalg.norm(ref))
        print(f"ns={r['ns']:4d} |J-J_ref|/|J_ref|={abs(r['J'] - rows[-1]['J']) / abs(rows[-1]['J']):.2e} "
              f"|g-g_ref|/|g_ref|={np.linalg.norm(g - ref) / np.linalg.norm(ref):.2e} cos={cos:.6f}")


if __name__ == "__main__":
    main(sys.argv[1:])
