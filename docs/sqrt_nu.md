# Strong-precession offset from tabulated coefficients

`sqrt_nu_rabe.x` evaluates the analytical WKB distribution by one-dimensional
quadrature for the SD2.2 strong-precession model. It is independent of the full
finite-precession boundary-value solver. The tabulated-coefficient interface is
intended for controlled model checks and later geometry integration. The companion `sqrt_nu_field.x` producer derives these tables from a
VMEC, boozmn or chartmap input for supported QI helicity. No finite-electric-field
DKES or MONKES benchmark is claimed.

Run `sqrt_nu_rabe.x input.in output.nc`, or omit the arguments to use
`sqrt_nu.in` and `sqrt_nu.nc`. The input contains two consecutive namelists:

```fortran
&sqrt_nu_sizes
 n_surface=2, n_eta=3, n_mode=1, n_nu=2, n_vE=3
/
&sqrt_nu_config
 s_tor=0.25,0.5
 nu_star=0.001,0.002
 v_E_star=-1.0,0.0,1.0
 eta(:,1)=1.0,1.5,2.0
 eta(:,2)=1.0,1.5,2.0
 tau_b_seconds=6*1.0
 tau_d_per_nu=6*1.0
 omega_theta_per_vE=6*-1.0
 dIv_deta=6*1.0
 g_modes=4*1.0
 normalization=2*1.0
 magnetic_drift_included=.false.
 eta_units='T-1'
 iv_units='synthetic example units'
/
```

This small synthetic example explains the storage order and is too coarse for
resolved WKB results. The registered `TestSqrtNuDriver` test uses a resolved
constant-coefficient model on two surfaces, two collisionalities and three
signed electric fields, and compares the result with an independently
integrated exponential expression. It also checks exact zero-precession
dispatch, output indexing and rejection of magnetic precession.

All geometry tables use the Fortran order `(eta,surface)`, except
`dIv_deta(eta,mode,surface)` and `g_modes(mode,collisionality,surface)`.
Harmonics are the positive integers `m=1,...,n_mode`. The caller supplies the
following quantities in a single consistent physical unit system:

| Input | Meaning and units |
| --- | --- |
| `eta` | Strictly increasing pitch grid from transition to deeply trapped endpoint, usually T^-1 |
| `tau_b_seconds` | Physical bounce time in seconds on that grid |
| `tau_d_per_nu` | `(tau_b*<D^etaeta>_b)/nu_star`, in eta squared |
| `omega_theta_per_vE` | Signed physical `Omega^theta/v_E_star`, in s^-1 |
| `dIv_deta` | Actual eta derivative of the sine coefficient `I_V^m` |
| `g_modes` | Boundary-layer harmonics evaluated at each requested collisionality |
| `normalization` | `pi*(dr/dpsi)/(integral dtheta0 integral dphi/B^2)` |

The diffusion table is multiplied by each `nu_star`, while the precession table
is multiplied by each signed `v_E_star`. The caller is responsible for the
linear collision scaling, the units of `g_modes` and `I_V`, and the physical
normalization that makes the final bootstrap offset dimensionless. The optional table `iv_deep(mode,surface)` supplies the actual `I_V^m` at
the final pitch point and defaults to zero. The kernel includes its boundary
term `gamma_s(deep)*iv_deep`, so a nonzero value must not be discarded. The
optional `iv_transition` and `iv_eta_b` tables record the first-grid-point
value and the native E77 boundary value, respectively.

The physical precession argument is `Omega^theta` in s^-1. The corresponding
library convention is `Omega_hat=Omega^theta/v`, in inverse length. With
`v_E_star=+c*E_r/(v*B_ref)`, the electric-drift expression gives
`Omega_hat=-v_E_star/r`, so a constant-radius input table has
`omega_theta_per_vE=-v/r`. The source note's equation labelled
`Omegatheta_cont` prints a plus sign in its equality to `v_E_star*v/r` while
also defining positive `v_E_star=c*E_r/(v*B_ref)`. Those two signs are
inconsistent. This interface follows the preceding explicit electric-drift
expression `Omega^theta=-c*E_r/(r*B_ref)` and records the convention in the
output. No undocumented sign reversal is applied to a supplied table.

Setting `magnetic_drift_included=.true.` stops before output is written.
Mixed zero and nonzero precession, or a precession sign change on a pitch grid,
is rejected by the kernel and recorded as a nonzero case status with NaN
results. The method does not resolve resonant or superbanana layers.

For exactly zero precession the kernel dispatches to the supplied E77
integral, rather than extrapolating the WKB approximation to its invalid
limit. This is an algebraic regression on the input coefficients, including the
supplied endpoint term. The native field tables start near `eta_c`, whereas
main RABE evaluates its E77 offset at `eta_b`. Their difference is retained
as a truncation diagnostic and must not be called exact recovery of the main
geometry calculation.

The NetCDF file has coordinate dimensions `surface`, `collisionality` and
`electric_field`. Each `case` maps to these coordinates through the one-based
`surface_index`, `nu_index` and `vE_index` variables. The `mode_sample`
dimension maps to a case and harmonic through `case_index` and `mode_index`.
It contains `flux_mode`, `lambda_mode`, `F_E_mode`, `g_mode`, `wkb_action` and
`valid_wkb`, together with `Iv_transition_mode` and `Iv_eta_b_mode`. Each
case also has `lambda_off`, `status`, `zero_precession` and the independently
assembled `lambda_E77_anchor`. The anchor is NaN when `iv_eta_b` is absent.
Units and definitions accompany the variables and global attributes.

`wkb_action` is the numerical integral of
`sqrt(abs(m*Omega^theta*tau_b/tau_d))` over eta. It is not the logarithmic
near-separatrix model parameter `A_m`, which requires additional local field
curvature information. The module exposes `compute_A_m` for that diagnostic
when the caller supplies the required geometry, and provides both the
logarithmic-model quadrature `fe_log_model` and closed-form
`fe_log_asymptotic` for comparison. The table driver does not infer these
geometry inputs from its numerical action. `valid_wkb` reports only the kernel's action and grid
resolution diagnostics. It does not establish near omnigeneity, mild rotation,
long mean free path or the validity of the boundary-layer harmonics. In
particular, no proposed domain threshold or finite-field benchmark tolerance
is promoted to a validated production acceptance criterion by this driver.


## Native QI field producer

Run the two stages in order:

```sh
sqrt_nu_field.x field.in tables.in
sqrt_nu_rabe.x tables.in result.nc
```

A field input has the same `sqrt_nu_sizes` group followed by this group:

```fortran
&sqrt_nu_field_config
 field_file='wout_published_QI.nc', field_type='vmec_nc'
 s_tor=0.25,0.5
 nu_star=0.0001,0.001
 v_E_star=-0.0001,0.0,0.0001
 radius_dkes=1.0,1.0
 speed=1.0
 M_pol=0.0, N_tor=4.0, max_n_fieldlines=32
 sign_sqrtg=-1, endpoint_cutoff=0.001
 magnetic_drift_included=.false.
/
```

Here `n_surface=2`, `n_nu=2` and `n_vE=3` must match the arrays. Set `n_eta`
and `n_mode` in the sizes group as well. The example radius values are
placeholders, not calibrated DKES radii for a particular equilibrium. Supply
physical `radius_dkes` values in metres and particle `speed` in m/s. The
producer uses `Omega^theta=-speed*v_E_star/radius_dkes`. The bootstrap gradient
normalization follows native RABE's effective radius `r_eff`, while the
precession conversion uses the supplied DKES radius. Their conversion must be
established before comparison with reference transport coefficients. A speed
of one is convenient because speed cancels between the bounce time and the
precession, but this does not determine the radius conversion.

Only `M_pol=0, N_tor=nfp` is supported, which makes RABE's sampled label
`xi_0` equal to the field-line label `theta_0` used by the precession equation.
Other helicities require a label transformation and are rejected. Each
neighbouring-maxima interval must contain exactly one minimum. Split maxima,
multiple minima and the existing native field-origin or symmetry sanity
failures are rejected. The producer does not disable those checks or repair
an input by an undocumented phase shift.

An explicitly calibrated reference well can instead be selected with
`reference_origin_theta` and `reference_origin_phi` (radians, both default zero)
in `sqrt_nu_field_config`. These are physical Boozer coordinates, not a change
to the field or input file. The producer places the label origin there and
checks the same stellarator symmetry, minimum and split-well requirements
about that reference. Its value must come from the input's angle convention,
not a search for passing WKB diagnostics. The two reference coordinates are
retained in the table and NetCDF output. An analytic cosine-well regression
checks the local pitch geometry, label placement and bounce/diffusion integrals
against an independent reference minimum while also retaining rejection of the
unmodified maximum at the zero origin.


The native `make_flock_of_fieldlines` calculation supplies the E77 boundary
integrals and the Fourier coefficients of `delta_aspect_ratio` and
`delta_eta`. With `l_c=pi*R/(2*nu_star)`, the generated assembly is

```text
g_m = <B^2>/<lambda_b> *
      [sqrt(eta_b*I_ref*l_c)*delta_aspect_m*S_A(m*iota_p)
       + (l_c/2)*delta_eta_m*S_B(m*iota_p)]
normalization = pi*(dr/dA_theta)/surface_average_normalization
```

Substitution into `normalization*sum(g_m*I_V^m(eta_b))` gives precisely
`Lambda_A/sqrt(nu_star)+Lambda_B/nu_star` from native `calc_offset_coefficients`.
The registered synthetic QI test checks this identity against that native
routine, independently of the new generated expression.

For each sampled pitch, turning points are found by bisection between the
single magnetic minimum and its neighbouring maxima. Existing substituted
quadrature integrates `local_radial_drift` to obtain `I_V(theta_0,eta)`, while
a generated integrand gives `tau_b=2*(B_phi+iota*B_theta)/speed *
integral dphi/[B*sqrt(1-eta*B)]`. Fourier transformation and second-order
finite differences give `dI_V^m/deta`. The producer averages the bounce time
over the actual perturbed field lines as an approximation to the note's
unperturbed coefficient. Its perturbative truncation error has not been
validated on general QI equilibria.

Following the boundary-frozen form of FEdef, the regular diffusion product is
held at `tau_d/nu_star=4*eta_b*mean(I_j(eta_b))/(pi*R)`, with
`I_j=(B_phi+iota*B_theta)*I_B` and the existing field-line `I_B` integrals.
This avoids assuming an unconfirmed eta dependence away from the boundary.
It remains a model choice, distinct from the numerical pitch discretization.

The pitch grid spans `eta_c+cutoff*(eta_t-eta_c)` to
`eta_t-cutoff*(eta_t-eta_c)`, where `eta_c=1/min(B_max)` and
`eta_t=1/max(B_min)` over the field lines. The explicit `endpoint_cutoff`
avoids the separatrix and deeply trapped endpoint singularities. An inward
turning-point guard of `1e-12` times the well's coordinate width prevents
roundoff from evaluating outside the trapped interval. Reduce the pitch
cutoff and refine both pitch and label grids when studying a configuration.
The synthetic QI test checks contraction on 33/65/129 pitch points and
16/32/64 field lines. It does not validate these defaults for a real field.

The supplied `iv_deep` preserves the finite-endpoint term. Comparing
`Iv_transition_mode` with `Iv_eta_b_mode` separates the eta-c versus eta-b
truncation from the WKB approximation. The independent `lambda_E77_anchor`
uses the original native eta-b integrals and is available at every scan
collisionality. Finite-field references and equilibrium-specific endpoint,
Fourier and WKB convergence studies are still required before production use.
