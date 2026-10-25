#!/usr/bin/env python3
"""F1 for NEO (SD2.2 criterion v3, section 2): sign of NEO lambda_b on a synthetic
tokamak-like Boozer file (M=1, N=0) against lambda_bB^tok = 1.46*sigma/(iota*sqrt(eps)).

NEO (github.com/itpplasma/NEO @ aea7a0c, calc_cur=1 with calc_eps=1) has no
sigma_sqrtg input. Measured (16 runs): sign(lambda_b) depends on sign(iota) and on the
handedness of the (R,Z) geometry (sign of Z), and not on sqrtg00 or on the sign of
curr_pol. With Z = +r sin(theta), R = R0 + r cos(theta) the Jacobian of
(r,theta,phi) -> (R,phi,Z) is negative, so sigma_sqrtg = -sign(Z). Hence sigma is
taken from the geometry here (this identification is a hypothesis; the code path is
isqrg = B**2/(G+iota*I) and psi_r, psi_z at neo_sub.f90:1873-1895).

usage: neo_f1_sign.py <neo.x executable> <work directory> [NEO source dir]
"""
import itertools
import math
import pathlib
import subprocess
import sys

NS, R0, A, B0 = 6, 5.0, 1.0, 2.0
S_INDEX = 3
NEO_SOURCE = pathlib.Path(
    sys.argv[3] if len(sys.argv) > 3 else "~/code/NEO"
).expanduser()


def write_bc(path, iota, pol_sign, sqrtg_sign, z_sign):
    rows = ["F1 synthetic tokamak", "h2", "h3", "h4", "h5", f"1 0 {NS} 1 1.0 {A} {R0}"]
    for i in range(NS):
        s = (i + 1) / NS * 0.5
        r = A * math.sqrt(s)
        eps = r / R0
        # file curr_pol -> G = -2e-7*curr_pol; pol_sign=+1 gives G = +R0*B0 < 0 flipped
        curr_pol = -pol_sign * R0 * B0 / 2e-7
        rows += [
            f"surf {i}",
            "hdr",
            f"{s:.8e} {iota:.8e} {curr_pol:.8e} 0.0 0.0 {sqrtg_sign * 1.0:.8e}",
            "m n r z l b",
            f"0 0 {R0:.8e} 0 0 {B0:.8e}",
            f"1 0 {r:.8e} {z_sign * r:.8e} 0 {-B0 * eps:.8e}",
        ]
    path.write_text("\n".join(rows) + "\n")


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


def run(exe, work, iota, pol_sign, sqrtg_sign, z_sign):
    work.mkdir(parents=True, exist_ok=True)
    for old in work.glob("t.*dat"):
        old.unlink()
    write_bc(work / "t.bc", iota, pol_sign, sqrtg_sign, z_sign)
    (work / "neo.in").write_text(neo_input(NEO_SOURCE))
    res = subprocess.run(
        [str(exe)], cwd=work, capture_output=True, text=True, timeout=600
    )
    if res.returncode != 0:
        raise SystemExit(f"NEO failed: {res.stdout[-400:]} {res.stderr[-400:]}")
    cols = (work / "t.cur.dat").read_text().split()
    return float(cols[2]), float(cols[5])  # lambda_b, lambda_b1


def main():
    exe, base = pathlib.Path(sys.argv[1]).resolve(), pathlib.Path(sys.argv[2]).resolve()
    eps = A * math.sqrt(S_INDEX / NS * 0.5) / R0
    print(
        "pol_sign iota_sign sqrtg00 Zsign | sigma=-sign(Z) | lambda_b "
        "lambda_b1 | sign(lambda_tok) s_N(lambda_b)"
    )
    failed = False
    baseline = {}
    s_n_all = set()
    for pol, iota_sign, sg, zs in itertools.product((1, -1), (1, -1), (1, -1), (1, -1)):
        iota = 0.5 * iota_sign
        lam_b, lam_b1 = run(
            exe, base / f"p{pol}_i{iota_sign}_g{sg}_z{zs}", iota, pol, sg, zs
        )
        # sigma_sqrtg = -sign(Z) (Jacobian of (r,theta,phi) -> (R,phi,Z)).
        sigma = -zs
        tok = 1.46 * sigma / (iota * math.sqrt(eps))
        s_n = math.copysign(1, lam_b) * math.copysign(1, tok)
        print(
            f"{pol:+d} {iota_sign:+d} {sg:+d} {zs:+d} | {sigma:+d} | {lam_b:+.6e} "
            f"{lam_b1:+.6e} | {math.copysign(1, tok):+.0f} {s_n:+.0f}"
        )
        baseline.setdefault((iota_sign, zs), set()).add(math.copysign(1, lam_b))
        s_n_all.add(s_n)
    for key, signs in baseline.items():
        if len(signs) != 1:
            print(f"sqrtg00 / curr_pol change sign(lambda_b) at iota,Z = {key}")
            failed = True
    print("s_N values over all runs:", sorted(s_n_all))
    if s_n_all != {1.0}:
        failed = True
    if failed:
        raise SystemExit(1)
    print("lambda_b sign is independent of sqrtg00 and curr_pol; s_N = +1")


if __name__ == "__main__":
    main()
