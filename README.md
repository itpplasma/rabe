# rabe

An implementation of the nea**r**-omnigenous, **a**symptotic **b**ootstrap **e**xpressions of Ref. [1].

## Summary

In case of the Lorentz collision model,
the mono-energetic bootstrap and Ware pinch coefficients, $\bar D_{31}$
and $\bar D_{13}$, can be expressed
via a dimensionless geometry parameter $\lambda_{bB}$,

$$
\bar D_{31}=\bar D_{13} = \frac{1}{3} v \rho_L B \lambda_{bB},
$$

scaling out particle velocity $v$, $\rho_L$ and $B$. This parameter can
be split into a ''symmetric'' contribution and a ''non-symmetric''  off-set,

$$
\lambda_{bB} = \lambda_\mathrm{sym} + \lambda_\mathrm{off}.
$$

For near-omnigenous stellarators at low plasma collisionality
$\nu_\ast=\pi R \nu_{\perp}/v$, where $\nu_\perp$ and $R$ are deflection frequency and
device major radius, respectively, the off-set can be expressed as
$$
    \lambda_{\text{off}} = \frac{\Lambda_\mathrm{A}}{\sqrt{\nu_\ast}} + \frac{\Lambda_\mathrm{B}}{\nu_\ast},
$$
where $\Lambda_\mathrm{A}$ and $\Lambda_\mathrm{B}$ are the geometrical factors
due to the variation of the trapped-passing boundary layer width and the misalignment
of local maxima, respectively (see Ref. [1]).

Given a VMEC equilibrium file, `rabe` outputs those geometric coefficients
needed to evaluate $\lambda_\mathrm{off}$ at any collisionality.

## Prerequisites

- Fortran compiler (e.g. gfortran)
- CMake ≥ 3.24
- [Ninja](https://ninja-build.org) build system
- [NetCDF-Fortran](https://github.com/Unidata/netcdf-fortran)
- BLAS and LAPACK
- [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse)

## Build and Tests

Use

```bash
make
```

to build the executable (debug build) or `make CONFIG=Release` for an optimized
build, and

```bash
make test
```

or

```bash
make test_all
```

to execute tests.

## Example

This example walks through a complete run using the QH stellarator equilibrium
from Ref. [2], which ships with the repository as a test input.

**Step 1 — build:**

```bash
make clean
make CONFIG=Release
```

This builds `rabe` in Release mode and writes the executable to `build/rabe`.

**Step 2 — create a working directory and copy the inputs:**

```bash
mkdir run_example && cd run_example
cp ../test/golden/input/rabe.in .
cp ../test/integration/vmec/input/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc .
```

For explanation of the input parameters see the Input / Output section below.

**Step 3 — run:**

```bash
../build/rabe
```

`rabe` reads `rabe.in` and `field_file` from the working directory and writes
results to `rabe.nc` after completion.

**Step 4 — visualize:**

```bash
cp ../test/golden/plot_golden.py .
python plot_golden.py rabe.nc
```

or with Octave:

```bash
cp ../test/golden/plot_golden.m .
octave plot_golden.m rabe.nc
```

Both produce `rabe_output.png` showing the off-set factors ($\Lambda_A$,
$\Lambda_B$), $\Lambda_\mathrm{HGM}$, and — if `should_calc_shaing_callen =
.true.` — the Shaing-Callen coefficients, all plotted against $s_\mathrm{tor}$.
The Python script requires `matplotlib` and `xarray`. The Octave script
requires the `netcdf` package (`pkg install -forge netcdf` if not present).

## Input / Output

`rabe` reads its configuration from a namelist file `rabe.in` in the working directory:

```fortran
&rabe_config
    field_file = "wout_example.nc",      ! VMEC equilibrium file (.nc)
    s_tor = 0.25, 0.5, 0.75,             ! surfaces by their normalized toroidal flux
    M_pol = -1.0,                        ! poloidal helicity of omnigenity
    N_tor = 4.0,                         ! toroidal helicity of omnigenity
    sign_sqrtg = -1.0,                   ! sign of the Jacobian sqrt(g)
    max_n_fieldlines = 200,              ! maximum field lines per surface
    should_calc_shaing_callen = .true.,  ! compute Shaing-Callen proxy
    n_eta = 100                          ! level resolution for Shaing-Callen
/
```

The type of omnigenity is set by the helicity of its perfect
omnigenous pendant, where contours of $B$ are closed after `M_pol`
toroidal and `N_tor` poloidal turns. `sign_sqrtg` is globally applied
to all coefficients to account for different coordinate conventions.
For typical `VMEC` output `sign_sqrtg=-1.0` (same as `signgs` in the `wou_*.nc`).

Results are written to `rabe.nc` (NetCDF), with one value per flux surface:

| Variable | Description |
|---|---|
| `s_tor` | normalized toroidal flux label |
| `Lambda_bl` | $1/\sqrt{\nu_\ast}$ factor ($\Lambda_\mathrm{A}$) |
| `Lambda_lm` | $1/\nu_\ast$ factor ($\Lambda_\mathrm{B}$)|
| `nu_star_crit` | lower collisionality limit for asymptotic model validity |
| `Lambda_finite` | $\sqrt{\nu_\ast}$ correction for finite boundary layer width |
| `err_flag` | 1 if omnigeneity violation is too strong, 0 otherwise |
| `R` | major radius [m] (reference length scale for $\nu_\ast$) |

`err_flag` warns the user to treat the results with caution,
as more than one local maximum on each side of a magnetic well was detected.

if `should_calc_shaing_callen = .true.`

| Variable | Description |
|---|---|
| `lambda_SC_bB` | omnigenous Shaing-Callen coefficient |
| `remainder` | non-omnigenous remainder of Shaing-Callen coefficient |

## Third Party

System libraries required at build time:

- [NetCDF-Fortran](https://github.com/Unidata/netcdf-fortran) for NetCDF output
- BLAS and LAPACK for linear algebra
- [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse) for sparse linear systems (needed for legacy `.bc` file field readers in tests)
- [libneo](https://github.com/itpplasma/libneo) for field file I/O (uses `$CODE/libneo` if present, otherwise fetched automatically)

Fetched automatically during build:

- [`quadpack`](https://github.com/jacobwilliams/quadpack) for numerical integration (BSD-3-Clause)
- [`pyplot-fortran`](https://github.com/jacobwilliams/pyplot-fortran) for visualization, source and license under `plot_lib` (BSD-3-Clause)

## References

[1] C.G Albert et al., *On the convergence of bootstrap current to the Shaing–Callen limit in stellarators*, Journal of Plasma Physics, 91(3), p. E77. [doi:10.1017/S0022377825000200](https://doi.org/10.1017/S0022377825000200) (2025)

[2] M. Landreman et al., *Optimization of quasi-symmetric stellarators with self-consistent bootstrap current and energetic particle confinement*, Phys. Plasmas 29, [doi:10.1063/5.0098166](https://doi.org/10.1063/5.0098166) (2022)

[3] A. Redl et al., *A new set of analytical formulae for the computation of the bootstrap current and the neoclassical conductivity in tokamaks*, Phys. Plasmas 28, [doi:10.1063/5.0012664](https://doi.org/10.1063/5.0012664) (2021)

[4] O. Sauter et al., *Neoclassical conductivity and bootstrap current formulas for general axisymmetric equilibria and arbitrary collisionality regime*, Phys. Plasma 6, [doi:10.1063/1.873240](https://doi.org/10.1063/1.873240) (1999)

[5] John R. Cary & Svetlana G. Shasharina, *Omnigenity and quasihelicity in helical plasma confinement systems, Phys. Plasmas 4*, 3323–3333, [doi:10.1063/1.872473](https://doi.org/10.1063/1.872473) (1997)
