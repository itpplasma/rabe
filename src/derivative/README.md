# Derivative mode: exact Jacobians with respect to the Boozer spectrum

The derivative mode computes, per flux surface, the offset coefficients
`Lambda_A`, `Lambda_B`, `Lambda_S`, `nu_star_crit` and `eta_b` together with
their exact Jacobian with respect to the Boozer spectrum of that surface. It is
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
- **Structure.** `find_structure` freezes the discrete choices with the
  standard algorithms: the continued-fraction approximation `p/q` of `iota_p`,
  the labels, and the pair of maxima that bound each line (bracket search and
  `pick_maximum_on_each_side`).
- **Field-line maxima.** They are polished by Newton on `g = B_phi + iota
  B_theta = 0` and not nudged inward. Their tangents follow from the implicit
  function theorem, `dphi* = -dg/(B_phiphi + 2 iota B_thetaphi + iota^2
  B_thetatheta)`, and `B_max` from the envelope theorem.
- **eta_b and nu_star_crit.** By default (`continuum_extrema = .true.`) the
  discrete maximum over labels is replaced by the global maximum of `B` on the
  surface, and the discrete minimum of the left maxima by the saddle of `B` on
  the maximum ridge (the mountain pass), both found by 2D Newton from the
  lattice argmax and argmin. These are the continuum limits of the discrete
  selections, their derivatives are the basis functions at the critical point
  (envelope theorem), and they do not depend on `iota`. The difference to the
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
   all spectrum directions.
4. **Moving endpoints and square-root singularities.** Lines whose maximum lies
   just below `1/eta_b` have near-singular derivatives of order `log(1/delta)`.
   They are integrable over the label, but a fixed rule resolves them only
   approximately, which is why `n_gauss = 64` is the default.
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
tangent and adjoint (`boozer_vmec_ad`, libneo PR
[#449](https://github.com/itpplasma/libneo/pull/449)). `bmnc`, `iota`, `G`
and `I` of each surface then come from libneo's transform, and
`rabe_deriv.nc` holds in addition the gradients of every output with respect to
the wout arrays `rmnc`, `zmns`, `lmns` (half-mesh, as stored in the wout) and
`iotaf`. For a cotangent `v[surface, output]` the wout cotangent is the
contraction `sum v * grad_X`, which is the VJP the VMEC++ adjoint needs; the
library routine `wout_chain_vjp` computes the same product with one libneo
adjoint sweep, and `wout_chain_jvp` the tangent.

- `nabla_s_mnc`, `psi_edge` and `R` still come from the projection of rabe's
  own field and are held fixed, so the gradients miss their dependence on the
  geometry. libneo's transform does not produce `|grad s|`.
- The wout arrays `bmnc`, `bsubumnc` and `bsubvmnc` do not enter libneo's
  transform, which rebuilds `B` from `R`, `Z`, `lambda` and `iota`, so their
  gradients are zero and are not written.
- The derivatives are those of libneo's discrete algorithm (splines in `rho`
  and the angles, Lagrange re-interpolation with a frozen stencil), not of the
  continuous Boozer transform; see the libneo PR for its kink at
  stellarator-symmetric points.
- The `iota` column stays the frozen-rational surrogate: chaining it with
  libneo's `iota` tangent moves the field-line iota with the true iota at fixed
  `p/q` and labels. libneo's `bmnc`, `G` and `I` also depend on `iota`, and
  that part is an ordinary derivative.

`TestWoutChainQH` and `TestWoutChainQA` perturb the wout harmonics in memory,
re-spline libneo's VMEC data and compare fourth-order central differences of
the composed primal (frozen structure, fixed `nabla_s_mnc`) with
`wout_chain_jvp`, separately for a geometry direction (`rmnc`, `zmns`, `lmns`)
and for `iotaf` at frozen field-line iota, and check the adjoint duality of
`wout_chain_jvp` and `wout_chain_vjp`.

## Usage

`rabe_deriv.x` reads `&rabe_config` from `rabe.in` (field file, helicity,
surfaces, `sign_sqrtg`, `max_n_fieldlines`) and the optional group

    &rabe_derivative
        mboz = 16, nboz = 16,        ! truncation of the libneo projection
        n_gauss = 64,                ! Gauss points per half well
        nu_star = 1e-3, 1e-2,        ! optional lambda_off evaluations
        continuum_extrema = .true.,
        wout_derivative = .false.,   ! chain to the wout harmonics (VMEC input)
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
| `values` | (surface, output) | `Lambda_A, Lambda_B, Lambda_S, nu_star_crit, eta_b` |
| `jacobian` | (surface, output, param) | `d values / d params` |
| `bmnc`, `nabla_s_mnc` | (surface, mode) | spectrum used |
| `xm`, `xn`, `nfp` | (mode), scalar | mode numbers, `xn` includes `nfp` |
| `iota`, `iota_input` | (surface) | rational field-line iota, input iota |
| `G`, `I`, `psi_edge`, `R` | (surface) | scalar parameters |
| `q`, `p`, `n_lines`, `n_split_lines` | (surface) | frozen structure |
| `diagnostics` | (surface, diagnostic) | gaps, Newton shift, `min lambda^2`, resonant weights, continuum shifts, iota kinks, mode flag |
| `spectral_residual`, `spectral_asymmetry` | (surface) | projection quality |
| `lambda_off`, `lambda_off_jacobian` | (surface, nu), (surface, nu, param) | when `nu_star` is given |
| `grad_rmnc`, `grad_zmns`, `grad_lmns` | (surface, output, radius, mn_mode) | `d values / d` wout array, with `wout_derivative` |
| `grad_iotaf` | (surface, output, radius) | `d values / d iotaf`, with `wout_derivative` |
| `wout_spectrum_difference` | (surface) | `max|bmnc_libneo - bmnc_projected|/B_00`, with `wout_derivative` |

## Tests

`TestSpectrumDerivative` checks central-difference step sweeps against the
Jacobian per parameter block and along random directions on a synthetic QH and
QA spectrum (best relative errors 1e-12 to 1e-9), the discrete mode, the
Gauss-Legendre moments and adjoint duality. `TestSpectrumDerivativePrimal`
compares the smooth primal with the standard pipeline on the same exact field
and checks the Gauss-Legendre convergence and the convergence of the continuum
shift. `TestSpectrumDerivativeQH` and `TestSpectrumDerivativeQA` repeat the
finite-difference check on spectra projected from the Landreman-Paul QH and QA
equilibria.
