"""Compare two RABE NetCDF output files for golden record testing."""

import sys

SKIP_CODE = 77

try:
    import numpy as np
    import xarray as xr
except ImportError:
    print("SKIP: xarray or numpy not installed")
    sys.exit(SKIP_CODE)


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <expected.nc> <current.nc>")
        sys.exit(1)

    expected_path = sys.argv[1]
    current_path = sys.argv[2]

    expected = xr.open_dataset(expected_path)
    current = xr.open_dataset(current_path)

    failed = False

    for var in expected.data_vars:
        if var not in current.data_vars:
            print(f"FAIL: variable '{var}' missing from current output")
            failed = True
            continue
        try:
            np.testing.assert_allclose(
                current[var].values,
                expected[var].values,
                rtol=1e-10,
                atol=1e-14,
            )
            print(f"PASS: {var}")
        except AssertionError as e:
            print(f"FAIL: {var}")
            print(f"      {e}")
            failed = True

    for var in current.data_vars:
        if var not in expected.data_vars:
            print(f"WARNING: variable '{var}' missing in golden record")

    expected.close()
    current.close()

    if failed:
        sys.exit(1)
    print("All variables match.")


if __name__ == "__main__":
    main()
