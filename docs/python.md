# Python interface

In addition to the Fortran executable, `rabe` ships a Python package with the
same name that exposes the core computation through a thin f90wrap binding.

## Installation

As a wheel (no compiler required if a matching binary is available):

```bash
pip install rabe
```

Or directly from source (requires the same Fortran prerequisites
listed in [building.md](building.md)):

```bash
pip install -e .
```

## Examples

Ready-to-run scripts are in the `python/` directory:

| Script | Field type | Use case |
| --- | --- | --- |
| `python/example.py` | `BoozerField` | Real VMEC equilibrium from a `.nc` file |
| `python/example_fourier.py` | `FourierField` | Analytical Fourier-mode field |

**`BoozerField`** (`boozer_field_t`): loads a VMEC NetCDF 3D equilibrium and
converts it to Boozer coordinates. This is the same field representation used
by the executable.

**`FourierField`** (`fourier_field_t`): builds a 2D field directly from
a flat list of Boozer Fourier modes plus surface values for covariant components.
Useful for benchmark and optimisation studies.

Both field types are accepted by the main downstream API via `FlockOfFieldlines`:

```python
from rabe.fieldline_mod import FlockOfFieldlines
flock = FlockOfFieldlines(max_n_fieldlines, iota, field, M_pol, N_tor, nfp)
lambda_a, lambda_b = flock.calc_offset_coefficients(R, dr_dAtheta)
nu_star_crit       = flock.calc_nu_star_crit(R)
```
