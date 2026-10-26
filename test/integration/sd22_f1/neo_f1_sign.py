#!/usr/bin/env python3
"""F1 for NEO (SD2.2 criterion v3, section 2): sign of NEO lambda_b on a synthetic
tokamak-like Boozer file (M=1, N=0) against lambda_bB^tok = 1.46*sigma/(iota*sqrt(eps)).

NEO (github.com/itpplasma/NEO @ aea7a0c, calc_cur=1 with calc_eps=1) has no
sigma_sqrtg input. sigma is therefore measured from the Boozer file itself: the sign
of the Jacobian det d(R, phi_cyl, Z)/d(s, theta_B, phi_B), evaluated with NEO's own
lab_swi=8 conventions (neo_sub.f90 @ aea7a0c: n -> n*nfp at :1362-1365, R = sum
rmnc cos(m theta - n phi), Z = sum zmns sin(...), phi_cyl = phi_B - (2 pi/nfp) l at
:1804-1830), with d/ds by central differences between neighbouring file surfaces.
The same function measures sigma for any lab_swi=8 file (option --bc).

A second part measures the magnitude: lambda_b1 / lambda_bB^tok over a scan in eps
and with B0 and R0 changed at fixed eps (reported, not asserted).

usage: neo_f1_sign.py <neo.x executable> <work directory> [NEO source dir]
       neo_f1_sign.py --bc <file.bc> [surface index (1-based)] ...
"""
import itertools
import math
import pathlib
import subprocess
import sys

NS, R0, A, B0 = 6, 5.0, 1.0, 2.0
S_INDEX = 3


def write_bc(path, iota, pol_sign, sqrtg_sign, z_sign, a=A, r0=R0, b0=B0):
    rows = ["F1 synthetic tokamak", "h2", "h3", "h4", "h5", f"1 0 {NS} 1 1.0 {a} {r0}"]
    for i in range(NS):
        s = (i + 1) / NS * 0.5
        r = a * math.sqrt(s)
        eps = r / r0
        curr_pol = -pol_sign * r0 * b0 / 2e-7
        rows += [
            f"surf {i}",
            "hdr",
            f"{s:.8e} {iota:.8e} {curr_pol:.8e} 0.0 0.0 {sqrtg_sign * 1.0:.8e}",
            "m n r z l b",
            f"0 0 {r0:.8e} 0 0 {b0:.8e}",
            f"1 0 {r:.8e} {z_sign * r:.8e} 0 {-b0 * eps:.8e}",
        ]
    path.write_text("\n".join(rows) + "\n")


def _is_int(tok):
    try:
        int(tok)
        return True
    except ValueError:
        return False


def read_bc(path):
    """Boozer file as read by NEO inp_swi=6/8: params line, then per surface a
    6-number data line followed by 6-column mode lines (m, n integer)."""
    nfp, surfaces = None, []
    for line in pathlib.Path(path).read_text().splitlines():
        tok = line.split()
        try:
            vals = [float(t) for t in tok]
        except ValueError:
            continue
        if nfp is None:
            if len(vals) == 7:
                nfp = int(vals[3])
            continue
        if len(vals) == 6 and _is_int(tok[0]) and _is_int(tok[1]):
            surfaces[-1]["modes"].append(vals)
        elif len(vals) == 6:
            surfaces.append(
                dict(
                    s=vals[0],
                    iota=vals[1],
                    curr_pol=vals[2],
                    curr_tor=vals[3],
                    modes=[],
                )
            )
    return nfp, surfaces


def _eval(surf, nfp, th, ph):
    """R, Z, p=phi_cyl-phi_B and their theta/phi derivatives (NEO lab_swi=8)."""
    out = [0.0] * 9  # R, Z, p, R_t, R_p, Z_t, Z_p, p_t, p_p
    for m, n, rc, zs, ls, _ in surf["modes"]:
        n = n * nfp
        arg = m * th - n * ph
        c, s = math.cos(arg), math.sin(arg)
        k = -2.0 * math.pi / nfp
        out[0] += rc * c
        out[1] += zs * s
        out[2] += k * ls * s
        out[3] += -m * rc * s
        out[4] += n * rc * s
        out[5] += m * zs * c
        out[6] += -n * zs * c
        out[7] += k * m * ls * c
        out[8] += -k * n * ls * c
    return out


def jacobian_sign(nfp, surfaces, k, nt=32, nph=16):
    """Measured sign of R*det d(R,phi_cyl,Z)/d(s,theta,phi) on surface k (0-based)."""
    lo, hi = surfaces[k - 1], surfaces[k + 1]
    ds = hi["s"] - lo["s"]
    signs, dets = set(), []
    for it in range(nt):
        for ip in range(nph):
            th = 2 * math.pi * it / nt
            ph = 2 * math.pi / nfp * ip / nph
            R, Z, p, Rt, Rp, Zt, Zp, pt, pp = _eval(surfaces[k], nfp, th, ph)
            a, b = _eval(lo, nfp, th, ph), _eval(hi, nfp, th, ph)
            Rs, Zs, ps = (b[0] - a[0]) / ds, (b[1] - a[1]) / ds, (b[2] - a[2]) / ds
            # rows: R, phi_cyl = phi + p, Z; columns: s, theta, phi
            det = (
                Rs * (pt * Zp - (1 + pp) * Zt)
                - Rt * (ps * Zp - (1 + pp) * Zs)
                + Rp * (ps * Zt - pt * Zs)
            )
            dets.append(R * det)
            signs.add(math.copysign(1, det))
    return signs, min(dets), max(dets)


def neo_input(neo_source):
    """doc/neo.in of the NEO checkout with the settings this check needs."""
    lines = (neo_source / "doc" / "neo.in").read_text().splitlines()
    subs = {
        "in_file": "t.bc",
        "no_fluxs": "1",
        "max_m_mode": "8",
        "max_n_mode": "2",
        "lab_swi": "8",
        "inp_swi": "8",
        "calc_eps": "1",
        "calc_cur": "1",
    }
    out = []
    for line in lines:
        words = line.split()
        key = words[1] if len(words) > 1 else ""
        if key in subs:
            line = f"{subs[key]:<20}{key}"
        elif words and words[0] == "16" and words[1] == "61":
            line = str(S_INDEX)
        out.append(line)
    return "\n".join(out) + "\n"


def run(exe, work, neo_source, **bc):
    work.mkdir(parents=True, exist_ok=True)
    for old in work.glob("t.*dat"):
        old.unlink()
    write_bc(work / "t.bc", **bc)
    (work / "neo.in").write_text(neo_input(neo_source))
    res = subprocess.run(
        [str(exe)], cwd=work, capture_output=True, text=True, timeout=600
    )
    if res.returncode != 0:
        raise SystemExit(f"NEO failed: {res.stdout[-400:]} {res.stderr[-400:]}")
    cols = (work / "t.cur.dat").read_text().split()
    nfp, surfaces = read_bc(work / "t.bc")
    sig, dmin, dmax = jacobian_sign(nfp, surfaces, S_INDEX - 1)
    if len(sig) != 1:
        raise SystemExit(f"Jacobian changes sign on the surface: {dmin} {dmax}")
    return float(cols[2]), float(cols[5]), sig.pop()


def sign_table(exe, base, neo_source):
    eps = A * math.sqrt(S_INDEX / NS * 0.5) / R0
    print("pol iota sqrtg00 Zsign | sigma_meas | lambda_b lambda_b1 | sgn(tok) s_N")
    failed, baseline, s_n_all = False, {}, set()
    for pol, isg, sg, zs in itertools.product((1, -1), (1, -1), (1, -1), (1, -1)):
        iota = 0.5 * isg
        lam_b, lam_b1, sigma = run(
            exe,
            base / f"p{pol}_i{isg}_g{sg}_z{zs}",
            neo_source,
            iota=iota,
            pol_sign=pol,
            sqrtg_sign=sg,
            z_sign=zs,
        )
        tok = 1.46 * sigma / (iota * math.sqrt(eps))
        s_n = math.copysign(1, lam_b) * math.copysign(1, tok)
        s_n1 = math.copysign(1, lam_b1) * math.copysign(1, tok)
        print(
            f"{pol:+d} {isg:+d} {sg:+d} {zs:+d} | {sigma:+.0f} | {lam_b:+.6e} "
            f"{lam_b1:+.6e} | {math.copysign(1, tok):+.0f} {s_n:+.0f}"
        )
        baseline.setdefault((isg, zs), set()).add(math.copysign(1, lam_b))
        s_n_all |= {s_n, s_n1}
    for key, signs in baseline.items():
        if len(signs) != 1:
            print(f"sqrtg00 / curr_pol change sign(lambda_b) at iota,Z = {key}")
            failed = True
    print("s_N values (lambda_b and lambda_b1) over all runs:", sorted(s_n_all))
    return failed or s_n_all != {1.0}


def magnitude_scan(exe, base, neo_source):
    print("a R0 B0 eps | lambda_b1 lambda_tok ratio | lambda_b ratio_b")
    cases = [
        (0.125, R0, B0),
        (0.25, R0, B0),
        (0.5, R0, B0),
        (1.0, R0, B0),
        (0.5, R0, 1.0),
        (1.0, 2 * R0, B0),
    ]
    for a, r0, b0 in cases:
        eps = a * math.sqrt(S_INDEX / NS * 0.5) / r0
        lam_b, lam_b1, sigma = run(
            exe,
            base / f"mag_a{a}_R{r0}_B{b0}",
            neo_source,
            iota=0.5,
            pol_sign=1,
            sqrtg_sign=1,
            z_sign=-1,
            a=a,
            r0=r0,
            b0=b0,
        )
        tok = 1.46 * sigma / (0.5 * math.sqrt(eps))
        print(
            f"{a} {r0} {b0} {eps:.5f} | {lam_b1:+.6e} {tok:+.6e} {lam_b1 / tok:.6f}"
            f" | {lam_b:+.6e} {lam_b / tok:.6f}"
        )


def main():
    if sys.argv[1] == "--bc":
        nfp, surfaces = read_bc(sys.argv[2])
        ks = [int(x) - 1 for x in sys.argv[3:]] or [len(surfaces) // 2]
        for k in ks:
            sig, dmin, dmax = jacobian_sign(nfp, surfaces, k)
            sf = surfaces[k]
            print(
                f"{sys.argv[2]} surface {k + 1} s={sf['s']:.6g} iota={sf['iota']:+.6g}"
                f" curr_pol={sf['curr_pol']:+.6g} sigma_meas={sorted(sig)}"
                f" R*det in [{dmin:.6g}, {dmax:.6g}]"
            )
        return
    exe, base = pathlib.Path(sys.argv[1]).resolve(), pathlib.Path(sys.argv[2]).resolve()
    neo_source = pathlib.Path(sys.argv[3] if len(sys.argv) > 3 else "~/code/NEO")
    neo_source = neo_source.expanduser()
    failed = sign_table(exe, base, neo_source)
    magnitude_scan(exe, base, neo_source)
    if failed:
        raise SystemExit(1)
    print("measured sigma: lambda_b sign independent of sqrtg00 and curr_pol; s_N = +1")


if __name__ == "__main__":
    main()
