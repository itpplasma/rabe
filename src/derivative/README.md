# Derivative mode: exact Jacobians with respect to the Boozer spectrum

The derivative mode computes, per flux surface, the offset coefficients
`Lambda_A`, `Lambda_B`, `Lambda_S`, `nu_star_crit` and `eta_b` and the
alpha-particle criteria `Lambda_alpha` and `Gamma_well` (see
[alpha_criteria.md](alpha_criteria.md)) together with their exact Jacobian
with respect to the Boozer spectrum of that surface. It is
opt-in and lives in its own library (`derivative_lib`) and executable
(`rabe_deriv.x`), so `rabe.x` and its outputs are unchanged.

## Parameters and outputs

A surface is described by `surface_spectrum_t`, and all differentiated inputs
are packed into one vector with `K` modes,

    params = [bmnc(1:K), nabla_s_mnc(1:K), iota, G, I, psi_edge, R]

where `B = sum_k bmnc(k) cos(xm(k) theta - xn(k) phi)` and likewise for
`|grad s|`, with `xn` including the number of field periods (booz_xform
convention). `G` and `I` are the covariant components `B_phi` and `B_theta` in
T m (booz_xform `Boozer_G`, `Boozer_I`), `psi_edge` is the toroidal edge flux
in Wb as used by rabe, and `R` is the major radius in m. The Jacobian has the
shape `jacobian(n_params, n_outputs)` with `jacobian(p, o) = d output_o / d
param_p`, so the JVP is `matmul(u, jacobian)` and the VJP with respect to the
spectrum is `matmul(jacobian, v)` (`spectrum_jvp`, `spectrum_vjp`). The
Jacobian is small and dense, and one forward sweep with all parameter
directions at once costs about one second per surface for `K = 545` on one core
(Release build).

The sign convention was checked against booz_xform for the Landreman-Paul QH
and QA equilibria at `s = 0.5`: the libneo projection agrees with booz_xform's
`bmnc_b` to 1.4e-5 relative with the same `(xm, xn)`, and `iota` and `G` agree
to 1e-6, so a booz_xform spectrum can be used directly.

## Method per component

The primal of the derivative mode follows `make_flock_of_fieldlines`,
`calc_offset_coefficients`, `calc_finite_boundary_layer_correction` and
`calc_nu_star_crit`, with the changes needed to make it a smooth function of
the parameters.

- **Field.** `spectral_surface_t` sums the Fourier series exactly, so `B` and
  `|grad s|` are linear in their coefficients and `dB/dbmnc(k)` is the basis
  function itself. There is no spline.
- **Origin.** rabe measures the field-line labels from a minimum of the ideal
  omnigenous field at `theta = phi = 0`. A quasi-isodynamic spectrum
  (`M_pol = 0`, `N_tor = nfp`) in booz_xform convention has its maximum there,
  and the bounding-maxima search stopped with "fewer than two maxima" (SQuID,
  Nautilus). When `suspect_omnigenous_origin_not_minimum` flags the origin,
  `find_structure` moves it by half a period (`phi_origin = pi/nfp`, or
  `theta_origin = pi` for `N_tor = 0`), the other stellarator-symmetric point,
  as the half-period shift of the QI `.bc` files of the bootstrap benchmark
  did. The field is evaluated at the shifted angles with the given
  coefficients, so the Jacobian stays with respect to the input spectrum.
- **Structure.** `find_structure` freezes the discrete choices with the
  standard algorithms: the continued-fraction approximation `p/q` of `iota_p`,
  the labels, and the pair of maxima that bound each line (bracket search and
  `pick_maximum_on_each_side`). Two maxima on one side of a line whose values
  agree within the symmetry violation are mirror images, and the one closest
  to the line origin bounds the well. The exact Fourier sum of a symmetric
  spectrum has zero violation, so the tolerance has a roundoff floor of 5e-13
  (relative). Without it the choice between two mirror maxima followed their
  last-bit difference, and a relative spectrum perturbation of 1e-14 moved
  `Lambda_S` by 1.4e-4 and `Lambda_A` by 1 % (LP QA, s = 0.5, booz_xform
  spectrum, two of 192 lines).
- **Field-line maxima.** They are polished by Newton on `g = B_phi + iota
  B_theta = 0` and not nudged inward. Their tangents follow from the implicit
  function theorem, `dphi* = -dg/(B_phiphi + 2 iota B_thetaphi + iota^2
  B_thetatheta)`, and `B_max` from the envelope theorem.
- **eta_b and nu_star_crit.** By default (`continuum_extrema = .true.`) the
  discrete maximum over labels is replaced by the global maximum of `B` on the
  surface, found by 2D Newton started from every line maximum within 1e-2
  (relative) of the lattice maximum, keeping the highest. A single start from
  the lattice argmax can stop at a lower local maximum that the lattice happens
  to sample at its peak (LP QA, s = 0.3: 1.1e-6 below the global maximum),
  which made the primal jump whenever the argmax line switched. Its
  derivative is the basis function at the maximum (envelope theorem) and does
  not depend on `iota`.
  The discrete minimum of the left maxima is replaced by the minimum over the
  continuous label `xi` of `F(xi)`, the bounding left maximum of the line `xi`
  as picked by the bracket search (`ridge_saddle`). `F' = (N_tor - iota
  M_pol)/nfp B_theta` by the envelope theorem, and its sign change from - to +
  is bracketed within one label spacing of every lattice local minimum within
  1e-2 of the lattice minimum and bisected to roundoff. The minimum is either
  a saddle of `B` on the maximum ridge (smooth minimum of one branch,
  derivative = basis function at the saddle) or a crossing of two maxima of
  equal height on one line (two branches, e.g. mirror images at a symmetry
  label, refined by Newton on `F_a = F_b`, derivative
  `(F_a' g_b - F_b' g_a)/(F_a' - F_b')` with `g` the derivative of a branch at
  fixed label, which includes `B_theta (phi - phi_0)` in the `iota`
  direction); `nu_crit_crossing` in `diagnostics` tells which. 2D Newton on
  `grad B = 0`, used before, converged to nearby maxima of `B` on the flat
  ridge of LP QA at s = 0.7, where the minimum is a crossing. If the minimum
  is not found (or `F` jumps because the bounding maximum disappears),
  `nu_star_crit` falls back to the discrete minimum with a warning, and the
  mode flag in `diagnostics` is 0. The difference to the
  discrete values shrinks like the square of the label spacing (synthetic QA:
  1.4e-5 relative in `eta_b` at 31 lines, 4.7e-7 at 171 lines). With
  `continuum_extrema = .false.` the discrete selections are used, and their
  derivative is the mean over the tie set, see below.
- **Line integrals.** The four integrals between the maxima use a fixed
  Gauss-Legendre rule (`n_gauss = 64` per half) in the substitution
  `phi = a + (b-a) u^2/2` on the left half and `phi = b - (b-a) u^2/2` on the
  right half, which regularises the square-root endpoints. The tangent is the
  exact derivative of this discrete rule, including the motion of the nodes
  with the endpoints and the dependence on `eta_b` through `lambda = sqrt(1 -
  eta_b B)`. On the synthetic cases `n_gauss = 64` and `128` agree to 1e-12.
- **Reductions.** `I_j`, `I_ref`, `delta_A`, `delta_eta`, the label DFT (with
  the `iota_p` phase shift of `delta_eta` differentiated as `dE_k = -i k E_k
  dshift`), the surface averages and the coefficient formulas are
  differentiated by hand in forward mode.

## Ill-posed parts and the surrogates used

These follow `rabe-derivatives/PHILOSOPHY.md` of the design notes.

1. **iota enters through a rational approximation.** rabe replaces `iota` by
   `p/q` and uses it as field-line iota, so the coded output is piecewise
   constant in the true `iota`. The `iota` column is the derivative with respect
   to the field-line iota at frozen `p/q` and frozen labels, i.e. the derivative
   of the result rabe would return if `iota` were exactly that rational.
   `rabe_deriv.x` therefore reports `iota` = `p/q` and keeps the input value as
   `iota_input`.
2. **S_A and S_B are discontinuous.** They are differentiated piecewise
   (`S_A' = 0.3`, `S_B' = 0`). At the rational `iota` the modes with `m iota_p =
   0 mod pi` sit on the jumps, the primal sets them to zero only when the
   floating-point remainder is below 1e-15, and the primal is discontinuous there
   by an amount proportional to their weight. Those modes get zero slope, and
   `resonant_weight_A` and `resonant_weight_B` report the share of `|rd_sin
   delta_cos|` they carry (synthetic QH: 3.9e-9). Finite differences in `iota`
   must therefore be taken off the rational point, as the tests do.
3. **Discrete max and min over labels.** At the rational `iota` the right
   maximum of one label and the left maximum of the label `iota_p` further are
   the same physical point, so the discrete max and min have ties whose
   `iota`-slopes differ, and the frozen-label surrogate has a kink there (LP QH,
   60 lines: four tied maxima with one-sided slopes of opposite sign). The
   continuum extrema remove this kink, which is why they are the default. In
   the discrete mode the tie set (within 1e-10) is averaged, which gives the
   central derivative, and `eta_b_iota_kink` and `nu_crit_iota_kink` report the
   spread. Ties between mirror images or identical points keep equal slopes in
   all spectrum directions. The continuum extrema are themselves a maximum
   (or minimum) over the local critical points of `B`, so they are kinked where
   two distinct local maxima or ridge minima exchange order. In near
   quasi-axisymmetry the ridge is almost flat and such exchanges are close: on
   LP QA at s = 0.3 the order of two local maxima 1.1e-6 apart changes within a
   wout perturbation of relative size 5e-6. The Jacobian is then a one-sided
   derivative. The same holds for the bounding maxima frozen in the
   structure: when two distinct (not mirror) maxima on one side of a line
   exchange order, the well and the primal jump, which is inherent in rabe's
   definition of the well. The frozen structure keeps the derivative on one
   side of such a jump, and `n_split_lines` counts the lines where it can
   happen.
4. **Moving endpoints and square-root singularities.** Lines whose maximum lies
   just below `1/eta_b` have near-singular derivatives of order `log(1/delta)`.
   They are integrable over the label, but a fixed rule resolves them only
   approximately, which is why `n_gauss = 64` is the default. On the line of
   the global maximum the nodes next to it can give `1 - eta_b B` below its
   exact lower bound `1 - c0 = 2 machine_eps` by roundoff of the Fourier sum
   (LP QH, s = 0.1, `n_gauss = 256`: -1.6e-15). Values within
   `64 machine_eps eta_b sum|bmnc|` below zero are raised to that bound, larger
   ones stop the evaluation.
5. **Derivative of the discretization.** The result is the exact derivative of
   the smooth discrete primal described above, not of rabe.x with QAG and nudged
   maxima. On the synthetic QH case both primals agree to 9e-8, and on LP QH
   (`max_n_fieldlines = 200`) rabe_deriv.x and rabe.x agree to 6e-4 in
   `Lambda_A`, 4e-4 in `Lambda_B`, 2e-6 in `Lambda_S` and 1e-3 in
   `nu_star_crit` at `mboz = nboz = 16`, with the difference shrinking with the
   truncation (to 1.4e-4, 5e-6, 9e-7, 3e-4 at 24 modes).
6. **Spectrum source.** Without `spectrum_file` the libneo spline field is
   projected onto the truncated spectrum; `spectral_residual` reports
   `max|B - B_spectrum|/B_00`. The projection itself is not differentiated.
   With `wout_derivative = .true.` the spectrum comes from libneo's
   differentiated transform instead, see below.
7. **lambda_off.** `lambda_off = Lambda_A/sqrt(nu_star) + Lambda_B/nu_star` is
   the above-critical asymptote and its gradient is only meaningful for `nu_star
   >= nu_star_crit`.

## Chain to the VMEC wout harmonics

With `wout_derivative = .true.` (VMEC input, no `spectrum_file`) the module
`wout_chain` composes the spectrum Jacobian with libneo's VMEC -> Boozer
tangent and adjoint (`boozer_vmec_ad`, libneo PRs
[#449](https://github.com/itpplasma/libneo/pull/449),
[#450](https://github.com/itpplasma/libneo/pull/450) and
[#451](https://github.com/itpplasma/libneo/pull/451)). Every spectrum
parameter of each surface then comes from libneo's transform: `bmnc`, the
`|grad s|` spectrum `nabla_s_mnc` (1/m), `iota`, `G`, `I`, `psi_edge` =
`phi_edge/(2 pi)` and `R` = `Rmajor_p`. `rabe_deriv.nc` holds in addition the
gradients of every output with respect to the wout arrays `rmnc`, `zmns`,
`lmns` (half-mesh, as stored in the wout) and `iotaf`, and with respect to the
scalars `phiedge` (a scaling of the whole `phi` profile) and `rmajor`
(`Rmajor_p`). For a cotangent `v[surface, output]` the wout cotangent is the
contraction `sum v * grad_X`, which is the VJP the VMEC++ adjoint needs; the
library routine `wout_chain_vjp` computes the same product with one libneo
adjoint sweep, and `wout_chain_jvp` the tangent.

- rabe's own field and projection are not used in this mode, so the wout is
  read once (LP QH, 5 surfaces, 200 lines, `mboz = nboz = 16`, Release:
  14.4 s to 5.7 s wall; the 25 wout VJPs take 0.2 s).
- `lambda_half_to_full` (`'linear'`, libneo's default, or `'cubic'`) selects
  libneo's half -> full mesh map of `lmns` for the VMEC field (both modes).
- `B`, `G`, `I` and `psi_edge` are proportional to the flux, so `Lambda_A`,
  `Lambda_B`, `Lambda_S` and `nu_star_crit` do not depend on `phiedge` and
  their `grad_phiedge` vanishes to roundoff; `eta_b` scales like `1/B`. `R`
  is read from the file, as rabe's field does, and depends only on `rmajor`.
- The wout arrays `bmnc`, `bsubumnc` and `bsubvmnc` do not enter libneo's
  transform, which rebuilds `B` from `R`, `Z`, `lambda` and `iota`, so their
  gradients are zero and are not written.
- The derivatives are those of libneo's discrete algorithm (splines in `rho`
  and the angles, Lagrange re-interpolation with a frozen stencil), not of the
  continuous Boozer transform; see the libneo PRs for its kink at
  stellarator-symmetric points.
- The `iota` column stays the frozen-rational surrogate: chaining it with
  libneo's `iota` tangent moves the field-line iota with the true iota at fixed
  `p/q` and labels. libneo's `bmnc`, `G` and `I` also depend on `iota`, and
  that part is an ordinary derivative.

`TestWoutChainQH` and `TestWoutChainQA` perturb the wout data in memory
(libneo's `boozer_ad_perturb`) and compare fourth-order central differences of
the composed primal (frozen structure, every spectrum parameter recomputed)
with `wout_chain_jvp`, separately for a geometry direction (`rmnc`, `zmns`,
`lmns`), for `iotaf` at frozen field-line iota, for `phiedge` (including the
invariance above) and for `rmajor`, check the adjoint duality of
`wout_chain_jvp` and `wout_chain_vjp`, and compare libneo's `nabla_s_mnc`,
`psi_edge` and `R` with rabe's own field.

## Usage

`rabe_deriv.x` reads `&rabe_config` from `rabe.in` (field file, helicity,
surfaces, `sign_sqrtg`, `max_n_fieldlines`) and the optional group

    &rabe_derivative
        mboz = 16, nboz = 16,        ! truncation of the libneo projection
        n_gauss = 64,                ! Gauss points per half well
        nu_star = 1e-3, 1e-2,        ! optional lambda_off evaluations
        continuum_extrema = .true.,
        n_eta = 32,                  ! trapped levels per well (Gamma_well)
        n_gauss_well = 32,           ! Gauss points per half bounce orbit
        wout_derivative = .false.,   ! chain to the wout harmonics (VMEC input)
        lambda_half_to_full = 'linear', ! libneo lambda map: 'linear' or 'cubic'
        spectrum_file = "",          ! read the spectrum instead of a field
        output_file = "rabe_deriv.nc"
    /

The spectrum file has the same variables as the output spectrum: dimensions
`surface` and `mode`, `nfp`, `xm(mode)`, `xn(mode)`, `s_tor(surface)`,
`bmnc(surface, mode)`, `nabla_s_mnc(surface, mode)` and `iota`, `G`, `I`,
`psi_edge`, `R` over `surface` (C order as seen from Python). An output file can
be read back as a spectrum file.

## Output layout (`rabe_deriv.nc`, dimensions in C order)

| variable | shape | content |
|---|---|---|
| `values` | (surface, output) | `Lambda_A, Lambda_B, Lambda_S, nu_star_crit, eta_b, Lambda_alpha, Gamma_well` |
| `jacobian` | (surface, output, param) | `d values / d params` |
| `bmnc`, `nabla_s_mnc` | (surface, mode) | spectrum used |
| `xm`, `xn`, `nfp` | (mode), scalar | mode numbers, `xn` includes `nfp` |
| `iota`, `iota_input` | (surface) | rational field-line iota, input iota |
| `G`, `I`, `psi_edge`, `R` | (surface) | scalar parameters |
| `q`, `p`, `n_lines`, `n_split_lines` | (surface) | frozen structure |
| `diagnostics` | (surface, diagnostic) | gaps, Newton shift, `min lambda^2`, resonant weights, continuum shifts, iota kinks, mode flag, `nu_crit_crossing` |
| `spectral_residual`, `spectral_asymmetry` | (surface) | projection quality (0 with `spectrum_file` or `wout_derivative`) |
| `lambda_off`, `lambda_off_jacobian` | (surface, nu), (surface, nu, param) | when `nu_star` is given |
| `grad_rmnc`, `grad_zmns`, `grad_lmns` | (surface, output, radius, mn_mode) | `d values / d` wout array, with `wout_derivative` |
| `grad_iotaf` | (surface, output, radius) | `d values / d iotaf`, with `wout_derivative` |
| `grad_phiedge`, `grad_rmajor` | (surface, output) | `d values / d phi_edge` (whole `phi` profile scaled), `d values / d Rmajor_p`, with `wout_derivative` |

## Tests

`TestSpectrumDerivative` checks central-difference step sweeps against the
Jacobian per parameter block and along random directions on a synthetic QH and
QA spectrum (best relative errors 1e-12 to 1e-9), the discrete mode, the
Gauss-Legendre moments and adjoint duality. `TestSpectrumDerivativePrimal`
compares the smooth primal with the standard pipeline on the same exact field
and checks the Gauss-Legendre convergence and the convergence of the continuum
shift. `TestSpectrumDerivativeQH` and `TestSpectrumDerivativeQA` repeat the
finite-difference check on spectra projected from the Landreman-Paul QH and QA
equilibria. `TestDerivativeRobustness_mirror` checks that 1e-14 spectrum noise
moves the outputs by less than 1e-9 on LP QA (s = 0.5, 37 lines with mirror
maxima; 1.2e-2 without the roundoff floor), `_saddle` that the ridge minimum
of LP QA at s = 0.7 is found and that `nu_star_crit` agrees with central
differences, and `_gauss256` that `n_gauss = 256` evaluates on LP QH at
s = 0.1 and agrees with 128 nodes.
