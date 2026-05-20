# rabe

An implementation of the nea**r**-omnigenous, **a**symptotic **b**ootstrap **e**xpressions of Ref. [1].

## Summary

In case of the Lorentz collision model,
the mono-energetic bootstrap and Ware pinch coefficients, $\bar D_{31}$
and $\bar D_{13}$, can be expressed
via a dimensionless geometry parameter $\lambda_{bB}$,

$$
\bar D_{31}=\bar D_{13} = \frac{1}{3} v \rho_L B \lambda_{bB}, \tag{1}
$$

scaling out particle velocity $v$, $\rho_L$ and $B$. This parameter can
be split into a ''symmetric'' contribution and a ''non-symmetric''  off-set,

$$
\lambda_{bB} = \lambda_\mathrm{sym} + \lambda_\mathrm{off}. \tag{2}
$$

For near-omnigenous stellarators at low plasma collisionality
$\nu_\ast=\pi R \nu_{\perp}/v$, where $\nu_\perp$ and $R$ are deflection frequency and
device major radius, respectively, the off-set can be expressed as
$$
    \lambda_{\text{off}} = \frac{\Lambda_\mathrm{A}}{\sqrt{\nu_\ast}} + \frac{\Lambda_\mathrm{B}}{\nu_\ast}, \tag{3}
$$
where $\Lambda_\mathrm{A}$ and $\Lambda_\mathrm{B}$ are the geometrical factors
due to the variation of the trapped-passing boundary layer width and the misalignment
of local maxima, respectively (see Ref. [1]).

Given a VMEC equilibrium file, `rabe` outputs those geometric coefficients
needed to evaluate $\lambda_\mathrm{off}$ at any collisionality.
Additionally, it can compute also the coefficients for the
''symmetric'' part

$$
\lambda_\mathrm{sym} = \lambda_{bB}^\mathrm{SC} + \lambda_{bB}^\mathrm{HGM}, \qquad \lambda_{bB}^\mathrm{HGM} = \Lambda_\mathrm{S}\sqrt{\nu_\ast} \tag{4}
$$

where the Shaing-Callen asymptotic $\lambda_{bB}^\mathrm{SC}$ is
approximated by the perfect omnigenous asymptotic
$\lambda_{bB}^\mathrm{SC} \rightarrow \lambda_{bB}^\mathrm{LC}$
obtained by Landreman and Catto [2] and
$\lambda_{bB}^\mathrm{HGM}$ is the contribution due to the finite
boundary correction derived by Helander, Geiger and Maasberg [3].

## Prerequisites

- Fortran compiler (e.g. gfortran)
- CMake >= 3.24
- [NetCDF-Fortran](https://github.com/Unidata/netcdf-fortran)
- BLAS and LAPACK
- [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse)

## Example

This example walks through a complete run using the QH stellarator equilibrium
from Ref. [4], which ships with the repository as a test input.
Commands are given for `bash` shell. While in `rabe`
run

**Step 1 — build:**

```bash
make clean
unset LIBNEO
make CONFIG=Release
```

This builds the executable `rabe.x` in Release mode and writes it to `build`.
The second line is only needed if you have an enviroment variable `LIBNEO` set
(see build instructions below).

**Step 2 — create a working directory and link the inputs:**

```bash
mkdir run_example && cd run_example
ln -s ../test/golden/input/rabe.in .
ln -s ../test/integration/vmec/input/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc .
```

For explanation of the input parameters see the Input / Output section below.

**Step 3 — run:**

```bash
../build/rabe.x
```

`rabe` reads `rabe.in` and `field_file` from the working directory and writes
results to `rabe.nc` (NetCDF) and `rabe.dat` (plain text) after completion.

**Step 4 — visualize:**

```bash
python ../test/golden/plot_golden.py rabe.nc
```

or with Octave:

```bash
octave ../test/golden/plot_golden.m rabe.nc
```

Both produce `rabe_output.png` showing the off-set factors ($\Lambda_A$,
$\Lambda_B$), $\Lambda_\mathrm{S}$, and — if `should_calc_shaing_callen =
.true.` — the Shaing-Callen asymptotic, all plotted against $s_\mathrm{tor}$.
The Python script requires `matplotlib` and `xarray`. The Octave script
requires the `netcdf` package (`pkg install -forge netcdf` if not present).

## Input / Output

`rabe` reads its configuration from a namelist file `rabe.in` in the working directory:

```fortran
&rabe_config
    field_file = "wout_example.nc",      ! VMEC equilibrium file (NetCDF)
    M_pol = -1.0,                        ! poloidal helicity of omnigenity
    N_tor = 4.0,                         ! toroidal helicity of omnigenity
    s_tor = 0.25, 0.5, 0.75,             ! explicit surface list OR
!   s_tor_min = 0.1,                     ! \
!   s_tor_max = 0.9,                     !  } n_s_tor equi-spaced surfaces
!   n_s_tor   = 9,                       ! /
    sign_sqrtg = -1.0,                   ! sign of the Jacobian sqrt(g)
    max_n_fieldlines = 200,              ! maximum field lines per surface
    should_calc_shaing_callen = .true.,  ! compute $\lambda_{bB}^\mathrm{LC}$
    n_eta = 100                          ! level resolution for trapped
                                         ! particle fraction computation
/
```

The type of omnigenity is set by the helicity

| type of omnigenity | `M_pol` | `N_tor` |
| --- | --- | --- |
| quasi-axisymmetric | 1 | 0 |
| quasi-isodynamic | 0 | $N_p$ |
| quasi-helicalsymmetric | $M$ | $N_p$ |

where $N_p$ is the number of field periods and $M$ is
given by the poloidal mode number $m$ of the strongest mode
$(m,n=N_p)$ of the Boozer spectrum

$$
 B(\vartheta, \varphi) = \sum_{m,n} B_{m,n} \cos{(m\vartheta - n\varphi)}. \tag{4}
$$

The list of surfaces on which to compute can either be given explicitly in `s_tor` OR via a uniform range. `sign_sqrtg` is globally applied
to all output coefficients to account for different coordinate conventions.
It should be set to the same value as `signgs` in the `VMEC` output `wout_*.nc` and is usually `sign_sqrtg=-1.0`.

Results are written to `rabe.nc` (NetCDF) and `rabe.dat` (plain text), with one
value per flux surface. Both files contain the same variables:

| Variable | Description |
| --- | --- |
| `s_tor` | toroidal flux $\psi$ normalized to flux at edge $\psi_a$ as $s_\mathrm{tor} = \psi/\psi_a$ |
| `Lambda_A` | $1/\sqrt{\nu_\ast}$ factor ($\Lambda_\mathrm{A}$ in Eq. 3) |
| `Lambda_B` | $1/\nu_\ast$ factor ($\Lambda_\mathrm{B}$ in Eq. 3) |
| `nu_star_crit` | lower collisionality limit for asymptotic model validity |
| `Lambda_S` | $\sqrt{\nu_\ast}$ correction for finite boundary layer width |
| `split_maxima` | 1 if omnigeneity violation is too strong, 0 otherwise |
| `R` | major radius [m] (reference length scale for $\nu_\ast$) |

`split_maxima` warns the user to treat the results with caution. If the
violation of omnigenity is too strong, local maxima contours are not merely
deformed, but also get split, which this flag notes.

if `should_calc_shaing_callen = .true.`, then

| Variable | Description |
| --- | --- |
| `lambda_LC_bB` | omnigenous Shaing-Callen coefficient $\lambda_{bB}^\mathrm{LC}$ |
| `remainder` | non-omnigenous remainder (only as a prototype) |

The `remainder` is a proxy of how much $\lambda_{bB}^\mathrm{LC}$ differs
from the actual $\lambda_{bB}^\mathrm{SC}$, but it does not include the
effect of bootstrap resonances and is not yet fully validated.

## Build and Tests

The build uses the enviroment variable `LIBNEO` to detect if a local copy of `libneo` (see third party dependencies below) is present. As a rule, if you are
not a developer, make sure to unset it before build e.g. with

```bash
unset LIBNEO
```

We use CMake for build configuration. If it is available on your machine, we
recommend [Ninja](https://ninja-build.org) as the generator

```bash
export CMAKE_GENERATOR=Ninja
```

Run

```bash
make
```

to build the executable (debug build) or `make CONFIG=Release` for an optimized
build. You can run the main suite of tests used to ensure correctness
of the code with

```bash
make test
```

or run all tests (including ones that take quite some time)

```bash
make test_all
```

The tests as well as their description, can be found in `test`.

## Third Party

System libraries required at build time:

- [NetCDF-Fortran](https://github.com/Unidata/netcdf-fortran) for NetCDF output
- BLAS and LAPACK for linear algebra
- [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse) for sparse linear systems (needed for legacy `.bc` file field readers in tests)

Fetched automatically during build:

- [libneo](https://github.com/itpplasma/libneo) for field file I/O; if a local copy at `$LIBNEO` is present, that one is used instead (MIT)
- [`quadpack`](https://github.com/jacobwilliams/quadpack) for numerical integration (BSD-3-Clause)
- [`pyplot-fortran`](https://github.com/jacobwilliams/pyplot-fortran) optional for visualization; source and license under `plot_lib` (BSD-3-Clause)

## References

[1] C.G Albert et al., *On the convergence of bootstrap current to the Shaing–Callen limit in stellarators*, Journal of Plasma Physics, 91(3), p. E77. [doi:10.1017/S0022377825000200](https://doi.org/10.1017/S0022377825000200) (2025)

[2] M. Landreman and P. J. Catto, *Omnigenity as generalized quasisymmetry*, Phys. Plasmas 19, 056103 [doi.org/10.1063/1.3693187](https://doi.org/10.1063/1.3693187) (2012)

[3] P. Helander, J. Geiger, and H. Maassberg, “On the bootstrap current in
stellarators and tokamaks”, Phys. Plasmas 18, 092505 [doi.org/10.1063/1.3633940](https://doi.org/10.1063/1.3633940) (2011)

[4] M. Landreman et al., *Optimization of quasi-symmetric stellarators with self-consistent bootstrap current and energetic particle confinement*, Phys. Plasmas 29, [doi:10.1063/5.0098166](https://doi.org/10.1063/5.0098166) (2022)



[5] John R. Cary & Svetlana G. Shasharina, *Omnigenity and quasihelicity in helical plasma confinement systems, Phys. Plasmas 4*, 3323–3333, [doi:10.1063/1.872473](https://doi.org/10.1063/1.872473) (1997)
