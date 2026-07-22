#!/usr/bin/env python3
"""Prove RABE's full output law under the pure VMEC relabel theta'=-theta."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

try:
    import netCDF4
    import numpy as np
except ImportError:
    print("SKIP: netCDF4 or numpy not installed")
    sys.exit(77)


SCALAR_SINE = {"rmns", "zmns", "bmns"}
ANGLE_COSINE = {"lmnc"}
ORIENTED_DENSITY_COSINE = {"gmnc"}
THETA_COMPONENT_COSINE = {"bsubumnc", "bsupumnc", "currumnc"}
UNCHANGED_COMPONENT_SINE = {"bsubvmns", "bsupvmns", "currvmns", "bsubsmns"}

SIGNED_OUTPUTS = (
    "Lambda_A",
    "Lambda_B",
    "Lambda_S",
    "lambda_LC_bB",
    "remainder",
)
SCALAR_OUTPUTS = ("s_tor", "nu_star_crit", "R")


def negate_if_present(dataset: netCDF4.Dataset, names: set[str]) -> None:
    for name in names:
        if name in dataset.variables:
            dataset.variables[name][:] = -dataset.variables[name][:]


def relabel_wout(source: Path, output: Path) -> None:
    """Write the same Cartesian equilibrium in the theta'=-theta chart."""

    shutil.copyfile(source, output)
    with netCDF4.Dataset(output, "a") as dataset:
        signgs = int(dataset.variables["signgs"][:])
        if signgs not in (-1, 1):
            raise RuntimeError(f"invalid source signgs={signgs}")
        dataset.variables["signgs"].assignValue(-signgs)
        for name in ("xn", "xn_nyq"):
            dataset.variables[name][:] = -dataset.variables[name][:]
        for name in ("iotaf", "iotas", "q_factor", "chi", "chipf"):
            if name in dataset.variables:
                dataset.variables[name][:] = -dataset.variables[name][:]

        negate_if_present(dataset, SCALAR_SINE)
        negate_if_present(dataset, ANGLE_COSINE)
        negate_if_present(dataset, ORIENTED_DENSITY_COSINE)
        negate_if_present(dataset, THETA_COMPONENT_COSINE)
        negate_if_present(dataset, UNCHANGED_COMPONENT_SINE)
        for name in ("buco", "jcuru"):
            if name in dataset.variables:
                dataset.variables[name][:] = -dataset.variables[name][:]

        dataset.setncattr(
            "rabe_manufactured_relabel",
            "theta_prime=-theta; zeta_prime=zeta; same Cartesian equilibrium",
        )


def input_text(field_file: str, m_pol: int, sign_sqrtg: int) -> str:
    return f"""&rabe_config
    field_file = \"{field_file}\",
    M_pol = {m_pol}.0d0,
    N_tor = 4.0d0,
    s_tor = 0.255102041, 0.316326531, 0.5, 0.704081633, 0.908163265,
    sign_sqrtg = {sign_sqrtg}.0d0,
    max_n_fieldlines = 200,
    should_calc_shaing_callen = .true.,
    n_eta = 100,
    unsafe_mode = .false.
/
"""


def run_rabe(executable: Path, directory: Path) -> subprocess.CompletedProcess[str]:
    environment = dict(os.environ)
    environment["OMP_NUM_THREADS"] = "1"
    return subprocess.run(
        [executable],
        cwd=directory,
        env=environment,
        check=False,
        text=True,
        capture_output=True,
    )


def relative_error(actual: np.ndarray, expected: np.ndarray) -> float:
    scale = max(float(np.linalg.norm(expected)), np.finfo(float).tiny)
    return float(np.linalg.norm(actual - expected) / scale)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rabe", type=Path, required=True)
    parser.add_argument("--wout", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        shutil.rmtree(output_dir)
    original_dir = output_dir / "original"
    relabel_dir = output_dir / "theta_relabel"
    mismatch_dir = output_dir / "wrong_explicit_sign"
    for directory in (original_dir, relabel_dir, mismatch_dir):
        directory.mkdir(parents=True)

    original_wout = original_dir / "wout.nc"
    relabel_wout_path = relabel_dir / "wout.nc"
    mismatch_wout = mismatch_dir / "wout.nc"
    shutil.copyfile(args.wout.resolve(), original_wout)
    relabel_wout(args.wout.resolve(), relabel_wout_path)
    shutil.copyfile(args.wout.resolve(), mismatch_wout)

    (original_dir / "rabe.in").write_text(input_text("wout.nc", -1, 0))
    (relabel_dir / "rabe.in").write_text(input_text("wout.nc", 1, 0))
    (mismatch_dir / "rabe.in").write_text(input_text("wout.nc", -1, 1))

    original_run = run_rabe(args.rabe.resolve(), original_dir)
    relabel_run = run_rabe(args.rabe.resolve(), relabel_dir)
    mismatch_run = run_rabe(args.rabe.resolve(), mismatch_dir)
    (original_dir / "run.log").write_text(original_run.stdout + original_run.stderr)
    (relabel_dir / "run.log").write_text(relabel_run.stdout + relabel_run.stderr)
    (mismatch_dir / "run.log").write_text(mismatch_run.stdout + mismatch_run.stderr)

    metrics: dict[str, float] = {}
    gates: dict[str, bool] = {
        "original_native_run_succeeds": original_run.returncode == 0,
        "relabel_native_run_succeeds": relabel_run.returncode == 0,
        "original_signgs_is_resolved": "resolved sign_sqrtg =           -1"
        in original_run.stdout,
        "relabel_signgs_is_resolved": "resolved sign_sqrtg =            1"
        in relabel_run.stdout,
        "explicit_wrong_sign_is_rejected": mismatch_run.returncode != 0
        and "inconsistent Boozer Jacobian sign"
        in (mismatch_run.stdout + mismatch_run.stderr),
    }
    if original_run.returncode == 0 and relabel_run.returncode == 0:
        with (
            netCDF4.Dataset(original_dir / "rabe.nc") as original,
            netCDF4.Dataset(relabel_dir / "rabe.nc") as relabeled,
        ):
            for name in SIGNED_OUTPUTS:
                error = relative_error(
                    np.asarray(relabeled[name][:]), -np.asarray(original[name][:])
                )
                metrics[f"{name}_opposite_relative_error"] = error
                gates[f"{name}_is_poloidal_covector"] = error < 5.0e-7
            for name in SCALAR_OUTPUTS:
                error = relative_error(
                    np.asarray(relabeled[name][:]), np.asarray(original[name][:])
                )
                metrics[f"{name}_same_relative_error"] = error
                gates[f"{name}_is_scalar"] = error < 5.0e-8
            gates["split_maxima_is_same"] = np.array_equal(
                original["split_maxima"][:], relabeled["split_maxima"][:]
            )

    evidence = {
        "coordinate_map": "theta_prime=-theta; zeta_prime=zeta",
        "physical_map": "same Cartesian position and magnetic field",
        "output_law": {
            "poloidal_covectors": list(SIGNED_OUTPUTS),
            "scalars": list(SCALAR_OUTPUTS) + ["split_maxima"],
        },
        "gates": gates,
        "metrics": metrics,
    }
    (output_dir / "vmec_theta_relabel_coefficients.json").write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(evidence, indent=2, sort_keys=True))
    if not all(gates.values()):
        raise RuntimeError(f"VMEC theta relabel coefficient gate failed: {gates}")


if __name__ == "__main__":
    main()
