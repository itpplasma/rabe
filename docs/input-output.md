# Input and output

`rabe` reads its configuration from a namelist file `rabe.in` in the working directory:

```fortran
&rabe_config
    field_file = "wout_example.nc",      ! VMEC equilibrium file (NetCDF)
!   field_type = "vmec_nc",              ! or "booz_xform", "chartmap"
    M_pol = -1.0,                        ! poloidal helicity of omnigenity
    N_tor = 4.0,                         ! toroidal helicity of omnigenity
    s_tor = 0.25, 0.5, 0.75,             ! explicit surface list OR
!   s_tor_min = 0.1,                     ! \
!   s_tor_max = 0.9,                     !  } n_s_tor equi-spaced surfaces
!   n_s_tor   = 9,                       ! /
    sign_sqrtg = -1.0,                   ! sign of the Jacobian sqrt(g)
    max_n_fieldlines = 200,              ! maximum field lines per surface
    should_calc_shaing_callen = .true.,  ! compute $\lambda_{bB}^\mathrm{LC}$
    n_eta = 100,                         ! level resolution for trapped
                                         ! particle fraction computation
    unsafe_mode = .false.                ! if .true., NaN-fill outputs for
                                         ! surfaces that fail a sanity check
                                         ! instead of aborting execution
/
```

`field_type` selects the format of `field_file`: a VMEC `wout` file
(`vmec_nc`, the default), a booz_xform output file (`booz_xform`) or a Boozer
chartmap (`chartmap`). `should_calc_shaing_callen` defaults to `.false.`.

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
 B(\vartheta, \varphi) = \sum_{m,n} B_{m,n} \cos{(m\vartheta - n\varphi)}.
$$

The list of surfaces on which to compute can either be given explicitly in `s_tor` OR via a uniform range. `sign_sqrtg` is globally applied
to all output coefficients to account for different coordinate conventions.
It should be set to the same value as `signgs` in the `VMEC` output `wout_*.nc` and is usually `sign_sqrtg=-1.0`.

By default (`unsafe_mode = .false.`), any failed sanity check halts the run
immediately with an error. This is the recommended behaviour as those checks point out not suited inputs e.g. violation of stellarator symmetry.
Setting `unsafe_mode = .true.` allows the run to continue, but outputs
are set to `NaN` for any surface where a check failed.

Results are written to `rabe.nc` (NetCDF) and `rabe.dat` (plain text), with one
value per flux surface. Both files contain the same variables:

| Variable | Description |
| --- | --- |
| `s_tor` | toroidal flux $\psi$ normalized to flux at edge $\psi_a$ as $s_\mathrm{tor} = \psi/\psi_a$ |
| `Lambda_A` | $1/\sqrt{\nu_\ast}$ factor ($\Lambda_\mathrm{A}$ in Eq. (3) of the [README](../README.md)) |
| `Lambda_B` | $1/\nu_\ast$ factor ($\Lambda_\mathrm{B}$ in Eq. (3)) |
| `nu_star_crit` | lower collisionality limit for asymptotic model validity |
| `Lambda_S` | $\sqrt{\nu_\ast}$ correction for finite boundary layer width |
| `split_maxima` | 1 if omnigeneity violation is too strong, 0 otherwise |
| `R` | major radius [m] (reference length scale for $\nu_\ast$) |

`split_maxima` warns the user to treat the results with caution. If the
violation of omnigenity is too strong, local maxima contours are not merely
deformed, but also get split, which this flag notes.

For fast runs, e.g. in optimization, we recommend to set
`should_calc_shaing_callen = .false.`. If it is enabled, the output
is extended by

| Variable | Description |
| --- | --- |
| `lambda_LC_bB` | omnigenous Shaing-Callen coefficient $\lambda_{bB}^\mathrm{LC}$ |
| `remainder` | non-omnigenous remainder (only as a prototype) |

The `remainder` is a proxy of how much $\lambda_{bB}^\mathrm{LC}$ differs
from the actual $\lambda_{bB}^\mathrm{SC}$, but it does not include the
effect of bootstrap resonances and is not yet fully validated.
