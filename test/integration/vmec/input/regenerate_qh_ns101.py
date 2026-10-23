"""Optional VMEC++ 0.9.1 recipe for the public QH ns101 fixture.

Usage: python regenerate_qh_ns101.py INPUT NEW_OUTPUT_DIRECTORY
Native tests use the tracked wout and never invoke this solver recipe.
"""
from pathlib import Path
import hashlib
import importlib.metadata as metadata
import json
import shutil
import sys
import time

import netCDF4
import numpy as np
import vmecpp

source, outdir = Path(sys.argv[1]), Path(sys.argv[2])
assert metadata.version("vmecpp") == "0.9.1"
assert hashlib.sha256(source.read_bytes()).hexdigest() == "1ccbdcdd29d79f13b02a49b3f139d628e556b036a41ea592236162f016327fde"
outdir.mkdir(parents=True, exist_ok=False)
work_input = outdir / "input.public_qh_ns101"
shutil.copyfile(source, work_input)
vi = vmecpp.VmecInput.from_file(str(work_input))
assert np.array_equal(vi.ns_array, [12, 101])
assert np.array_equal(vi.ftol_array, [1e-17, 1e-18])
assert np.array_equal(vi.niter_array, [1200, 20000])
settings = {}
for name in ("pmass_type", "am", "pcurr_type", "ac", "pres_scale", "gamma", "ncurr", "curtor", "phiedge", "nfp", "mpol", "ntor", "lasym", "lfreeb", "ns_array", "ftol_array", "niter_array"):
    value = getattr(vi, name)
    settings[name] = value.tolist() if hasattr(value, "tolist") else value
settings["unsupported_selector_attributes_present"] = {name: hasattr(vi, name) for name in ("pt_type", "ph_type", "at", "ah")}
(outdir / "converted-settings.json").write_text(json.dumps(settings, indent=2) + "\n")
vi.return_outputs_even_if_not_converged = True
start = time.monotonic()
output = vmecpp.run(vi, max_threads=1, verbose=False)
wout = outdir / "wout.nc"
output.wout.save(str(wout))
with netCDF4.Dataset(wout) as ds:
    result = {key: float(ds[key][...]) for key in ("ns", "ier_flag", "fsqr", "fsqz", "fsql", "ftolv")}
result.update(wall_s=time.monotonic()-start, vmecpp=metadata.version("vmecpp"), wout_sha256=hashlib.sha256(wout.read_bytes()).hexdigest(), wout_bytes=wout.stat().st_size)
result["converged"] = result["ns"] == 101 and result["ier_flag"] == 0 and all(np.isfinite(result[k]) and 0 <= result[k] <= 1e-17 for k in ("fsqr", "fsqz", "fsql"))
(outdir / "result.json").write_text(json.dumps(result, indent=2)+"\n")
print(json.dumps(result))
if not result["converged"]:
    raise SystemExit("QH solve did not meet the frozen convergence contract")
