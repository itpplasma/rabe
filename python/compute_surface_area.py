"""Compute flux surface areas from a VMEC wout file using simsopt.

Prints Fortran-formatted reference values for use in integration tests.
"""

import numpy as np
from simsopt.geo import SurfaceRZFourier

wout_file = (
    "test/integration/vmec/input/"
    "wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"
)

s_values = [0.1, 0.3, 0.5, 0.7, 0.9]

print(f"{'s':>5s}  {'area [m^2]':>24s}")
print("-" * 32)

areas = []
for s in s_values:
    surf = SurfaceRZFourier.from_wout(wout_file, s=s, ntheta=128, nphi=128)
    area = surf.area()
    areas.append(area)
    print(f"{s:5.1f}  {area:24.16e}")

print()
print("Fortran array (copy into test):")
print("    area_ref = [", end="")
for i, a in enumerate(areas):
    sep = ", &" if i < len(areas) - 1 else "]"
    nl = "\n                " if i < len(areas) - 1 else ""
    print(f"{a:.16e}_dp{sep}{nl}", end="")
print()
