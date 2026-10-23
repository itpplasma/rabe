# Alpha-particle criteria of the derivative mode: `Lambda_alpha` and `Gamma_well`

The derivative mode returns four outputs in addition to the offset
coefficients, `Lambda_alpha` (index 6), `Gamma_well` (index 7), `var_dJ`
(index 8) and `Gamma_s` (index 9), with exact gradients like the others. They implement two of the smooth optimization criteria that S. V.
Kasilov proposed for alpha confinement in analogy to the bootstrap offset.
Omnigeneity satisfies all of them:

1. align the maxima of `B` (`nu_star_crit` measures their spread),
2. at the trapped-passing boundary `eta = eta_b` make the integrals between
   the bounding maxima independent of the field-line label,
   (a) `int dphi lambda/B^phi` (no bootstrap offset, `Lambda_A`) and
   (b) `int dphi B lambda/B^phi` (no alpha losses, `Lambda_alpha`),
3. keep the bounce-averaged drift within the wells small, so that the
   superbanana orbits stay radially bounded (`Gamma_well`).

Criterion 2b is a proxy built from the offset machinery, and criterion 3 is a
proxy of the Gamma_c type. Neither is derived from alpha-particle transport
theory, and their correlation with alpha losses has to be checked with an orbit
code such as SIMPLE. Prompt losses of deeply trapped particles are outside both.

## Notation

Boozer angles `(theta, phi)`, field-line label `alpha = theta - iota phi`, line
`theta = theta_0 + iota (phi - phi_0)`, and `lambda = sqrt(1 - eta B) =
|v_par|/v` with the pitch variable `eta = v_perp^2/(v^2 B)`. In Boozer
coordinates with rabe's flux normalization `B^phi = B^2/(G + iota I)`, so that
the arc length is `dl = B dphi/B^phi = (G + iota I) dphi/B`.

## (2b) Lambda_alpha

The integral of criterion 2a at `eta_b` is rabe's

    I_j = int_{phi_1}^{phi_2} dphi lambda_b/B^phi = (G + iota I) int dphi lambda_b/B^2,

`lambda_b = sqrt(1 - eta_b B)`, between the bounding maxima `phi_1, phi_2` of
line `j` (Eq. 26 of Albert et al., J. Plasma Phys. 91, 2025, here `I_j`). It
sets the boundary-layer width `delta eta_j ~ I_j^(1/2)` and the aspect ratio
`A_j = delta eta_ref/delta eta_j` of the offset domains (Eqs. 84 and 88 there).
rabe uses `delta A_j = sqrt(I_ref/I_j) - 1` with `I_ref = (n/sum_j
I_j^(-1/2))^2`, which makes the label mean of `sqrt(I_ref/I_j)` one, and

    Lambda_A = pi sum_m rd_s(m) dA_c(m) S_A(m iota_p) * <B^2>/<lambda_b>
               * sqrt(eta_b I_ref)/norm * dr/dA_theta * sqrt(pi R/2),

where `rd_s` are the sine coefficients over the label of the radial-drift
integral, `dA_c` the cosine coefficients of `delta A_j`, and `S_A` the fit of
Eq. 133 (slope 0.3).

The integral of criterion 2b is the parallel adiabatic invariant of the
transitioning (barely trapped) particles,

    J_j = int_{phi_1}^{phi_2} dphi B lambda_b/B^phi = (G + iota I) int dphi lambda_b/B
        = int dl v_par/v  at eta = eta_b.

Criterion 2b asks for `dJ_j/dtheta_0 = 0`, so that the transitioning particles
see constant `J` on the surface and their bounce-averaged radial drift, which
is proportional to `dJ/dalpha` (Velasco et al., Nucl. Fusion 61, 116059, 2021,
Eq. 14), vanishes. `Lambda_alpha` treats the label variation of `J_j` exactly
like that of `I_j`,

    delta J_j = sqrt(J_ref/J_j) - 1,   J_ref = (n/sum_j J_j^(-1/2))^2,
    Lambda_alpha = Lambda_A with dA_c replaced by the cosine coefficients of delta J_j,

with the same prefactor (including `sqrt(I_ref)`), so that `Lambda_alpha` and
`Lambda_A` are directly comparable. The factor `G + iota I` cancels in
`delta J_j`. Two properties follow from the definition. First, `Lambda_alpha`
is a projection: it vanishes when `J_j` does not depend on the label, but it
can also vanish for a varying `J_j` whose cosine modes are orthogonal to `rd_s
S_A`. In an objective it is therefore used squared and together with the
other criteria. Second, with `B = B_0 (1 + O(eps))` along the well `J_j =
B_0 I_j (1 + O(eps))`, so close to omnigeneity `delta J_j` and `delta A_j`
share their leading label dependence and `Lambda_alpha = Lambda_A (1 +
O(eps))` (LP QH, `s = 0.3`: `4.03e-7` against `3.77e-7`).

## (3) Gamma_well

For a level `eta` trapped in the well of line `j`, between the turning points
`phi_-`, `phi_+` with `1 - eta B = 0`,

    J(eta, alpha) = (G + iota I) int_{phi_-}^{phi_+} dphi lambda/B,
    dJ/dalpha     = (G + iota I) int dphi d/dtheta (lambda/B)
                  = -(G + iota I) int dphi B_theta (1 + lambda^2)/(2 lambda B^2).

The end points contribute nothing because `lambda/B` vanishes there. The
bounce-averaged radial drift is proportional to `dJ/dalpha` (Velasco et al.
2021, Eq. 14), and Nemov et al., Phys. Plasmas 15, 052501 (2008) define
`Gamma_c` from `gamma_c = (2/pi) arctan(v_r/v_theta)`, the ratio of the radial
and the poloidal bounce-averaged drift, which in Velasco's notation is `(2/pi)
arctan(dJ/dalpha / |dJ/ds|)`. The radial derivative `dJ/ds` needs a second
surface, which the per-surface derivative mode does not have, so
`Gamma_well` normalizes `dJ/dalpha` by `J` instead,

    Gamma_well = sum_{j,k} W_jk (dJ/dalpha)_jk^2 / sum_{j,k} W_jk J_jk^2,

over the wells `j` of the flock (one per line, between its bounding maxima,
which covers the surface) and their trapped levels

    eta_jk = eta_lo,j + (eta_hi,j - eta_lo,j) t_k,   W_jk = (eta_hi,j - eta_lo,j) w_k,

with `eta_lo = 1/min(B_max,1, B_max,2)` (barely trapped), `eta_hi = 1/B_min`
(bottom), and `t_k = u_k^2`, `w_k = 2 u_k w_k^GL` from Gauss-Legendre on
`(0, 1)`. `Gamma_well` is dimensionless, invariant under `B -> c B` (so its
`phiedge` gradient vanishes, as for `Lambda_A`), and
vanishes for quasi-symmetric and omnigenous fields. A synthetic quasi-helical
field `B = 1 - 0.1 cos(theta + 4 phi) + eps cos(theta)` gives `Gamma_well =
1e-25` at `eps = 0` and grows like `eps^2` (pyrabe `test_alpha_criteria.py`).
Unlike `gamma_c` the integrand is not bounded: `dJ/dalpha` diverges like `log
t` at the barely trapped end, where the bounce time diverges, which the
nodes `t = u^2` integrate without loss.

## (2b, direct) var_dJ

`Lambda_alpha` is the projection of `delta J_j` onto the bootstrap weights
`rd_s S_A` and therefore scales like `eps^2` and can vanish by cancellation
of modes. The direct measure of criterion 2b is the variance of `delta J_j`
over the labels of the flock (output 8),

    var_dJ = (1/n) sum_j delta J_j^2 = sum over label modes of the squared
             cosine and sine coefficients of delta J_j (Parseval),

where the label mean of `delta J_j` vanishes by the choice of `J_ref`. Since
`delta J_j = -(J_j - <J>)/(2 <J>) + O(delta^2)`, `var_dJ` is a quarter of the
relative variance of `J` of the transitioning particles over the field-line
label. It vanishes only if `J_j` is the same on all lines and grows about
like `eps^2` near omnigeneity (`delta J_j ~ eps`; local exponent 1.93 between
`eps = 1e-4` and `1e-3` and 1.89 between `1e-3` and `1e-2` on the helical
family of `Gamma_well`). It is computed directly from
`delta J_j` and its tangent, which equals the Parseval sum of the label
transform exactly.

## (3, Gamma_c-like) Gamma_s

Criterion 3 asks that the superbanana orbits of the trapped particles stay
radially bounded, not that `J` is independent of the label inside the wells,
which `Gamma_well` (normalized by `J`) enforces. In the coordinates `(s,
alpha)` the bounce-averaged orbit is a contour of `J(s, alpha)` at fixed
`eta` and energy, `ds/dt ~ dJ/dalpha` and `dalpha/dt ~ -dJ/ds`, so the
radial excursion over a variation `Delta_alpha J` of `J` along the surface is
`Delta s ~ Delta_alpha J/|dJ/ds|` (Velasco et al. 2021, Eq. 14, for the
drifts). The criterion is therefore the ratio
`(dJ/dalpha)/(dJ/ds)`, which is `tan` of the angle between the drift and the
surface and is the quantity whose `arctan` Nemov et al. (2008) average in
`Gamma_c`. Output 9 is the mean-square form of that ratio,

    Gamma_s = sum_{j,k} W_jk (dJ/dalpha)_jk^2
              / sum_{j,k} W_jk ((dJ/dalpha)_jk^2 + (dJ/ds)_jk^2),

with the wells, levels and weights `W_jk` of `Gamma_well`. It is `Gamma_well`
with the normalization `J^2` replaced by `|grad_(s,alpha) J|^2`, lies in
`[0, 1]`, is invariant when `B`, its radial derivative, `G` and `I` are
scaled together, vanishes for quasi-symmetric and omnigenous fields, and grows
like `eps^2` (local exponent 1.98 between `eps = 1e-3` and `3e-3` and 1.93 up
to `1e-2` on the helical family). For small ratios it is the weighted mean of `(dJ/dalpha)^2/(dJ/ds)^2` with the
weights `W (dJ/ds)^2`.

The per-level form `sum W gamma/sum W`, `gamma = (dJ/dalpha)^2/((dJ/dalpha)^2
+ (dJ/ds)^2)`, which is closer to `Gamma_c` (and the plain ratio, which
diverges), was tried first and rejected. `dJ/ds` changes sign at the
precession reversal (maximum-`J` levels), where `gamma` has a peak of height
one and width `|dJ/dalpha|/|d^2J/(ds deta)|` in `eta`. Near symmetry that
width is far below the spacing of the trapped levels, so the sum did not
converge in `n_eta` (LP QH `s = 0.3` with a synthetic radial block: `1.3e-3`,
`2.3e-3`, `8.5e-3`, `2.7e-3` for `n_eta = 16, 32, 64, 128`) and was neither
smooth nor of a definite order in `eps` (`~eps^1.05`), whereas the ratio of
the sums converges (`2.921332e-5` for `n_eta = 32, 64, 128`). The price is
that a reversal level enters only with its weight, not with `gamma = 1`.

The radial derivative at fixed `alpha`, `eta` and energy is, with `G_t = G +
iota I` and primes for `d/ds`,

    dJ/ds = G_t' j + G_t j_s,   j = int dphi lambda/B,
    j_s   = -int dphi (1 + lambda^2)/(2 lambda B^2) (B_s + iota' (phi - phi_min) B_theta),

where `B_s = sum_mn dbmnc/ds cos(m theta - n phi)` is the radial derivative
of the Boozer spectrum at fixed Boozer angles, `G_t' = G' + iota' I + iota
I'`, and the end points contribute nothing because `lambda/B` vanishes there.
The term with `iota'` is the magnetic shear: at fixed label the line moves
in `theta` by `iota' (phi - phi_c) ds`, and `alpha` is measured at `phi_c =
phi_min`, the bottom of the well (its tangent follows from the implicit
function theorem on `B_phi + iota B_theta = 0`). A different `phi_c` adds
`(phi_min - phi_c) iota' dJ/dalpha` to `dJ/ds`, which changes `Gamma_s` by a
relative `O(eps)`. Since `G_t` cancels in `dJ/dalpha` but not in
`dJ/ds`, the sums use `D = dJ/dalpha/G_t` and `P = dJ/ds/G_t = j_s + (G_t'/G_t)
j`.

The radial derivatives are parameters of the surface: `params` gains the
block `dbmnc_ds(1:K), diota_ds, dG_ds, dI_ds` after `R` (`n_params = 3K + 8`),
on which only `Gamma_s` depends. They come from the source of the spectrum:

- `wout_derivative` (libneo's differentiated VMEC -> Boozer transform) and
  pyrabe's `WoutChain`: central differences of the Boozer spectra at
  `s +- ds_radial` (`ds_radial = 0.01` by default), each computed by the same
  transform. The Boozer angles are unique for a stellarator-symmetric
  surface, so the difference is `dB/ds` at fixed Boozer angles up to
  `O(ds^2)` (pyrabe, LP QH, `s = 0.25`: `Gamma_s` changes by `2.2e-8`,
  `5.4e-9`, `1.3e-9` from `ds_radial = 0.02` to `0.0025`, relative `3e-4` at
  the default). The gradient with respect to the wout is the exact gradient of
  this three-surface definition (two more transforms and their adjoints per
  surface). libneo's transform interpolates the half-mesh harmonics radially
  at a given `s`, but it is not differentiated with respect to `s` itself,
  which is why the radial derivative is a difference of neighbouring surfaces.
- rabe's own field (`spectralize_field`): the same central difference of the
  spectra of `s +- ds_radial`.
- `spectrum_file`: the optional variables `dbmnc_ds`, `diota_ds`, `dG_ds`,
  `dI_ds`. Without them `Gamma_s` is set to zero with zero gradient and a
  warning is printed.

## Numerical scheme

- Well minimum: lowest of 1025 samples between the polished bounding maxima,
  polished by Newton on `B_phi + iota B_theta = 0`; `dB_min` by the envelope
  theorem.
- Turning points: first crossing of `B = 1/eta` going outward from the
  minimum, bracketed on the samples and converged by safeguarded Newton
  (Fortran) or bisection plus Newton (pyrabe); the tangent from `1 - eta
  B(phi*) = 0`, `dphi* = -(B deta + eta dB|_phi)/(eta B')`. With a secondary
  minimum this is the turning point of the particle at the bottom.
- Integrals: `n_gauss_well` Gauss points per half in `phi = a + (b - a)
  u^2/2` from each end, which makes both integrands regular at the simple
  zeros of `lambda^2`. At a node `a + d`, `lambda^2 = -eta (B(a + d) - B(a))`
  with the cosine differences summed as `-2 sin(arg_a + c) sin(c)`, `c = (m
  iota - n) d/2`. The naive `1 - eta B` loses the digits of `B` next to the
  end, which made the primal noisy at the relative level `eps/u_1^2` (finite
  differences of `Gamma_well` stalled at `1e-5`; with the difference form they
  reach `8e-9`, LP QA, `mboz = nboz = 8`, iota direction).
- Defaults `n_eta = 32`, `n_gauss_well = 32`: on the LP QH and QA fixtures
  `n_gauss_well = 32` and `64` agree to `1e-10`, and `n_eta = 32` and `64` to
  `3e-7` (QH) and `2e-4` (QA), while `n_eta = 16` was 5 % off on QA.

The algebraic identities used above (the `theta` derivative of `lambda/B`, the
tangent of `(1 + lambda^2)/(2 lambda)`, the cosine difference, the tangent of
`delta J_j` and of the turning point) are proved by SymEngine through fortsym
in `symbolic/check_alpha_criteria.f90`.

## Kinks

`Gamma_well` inherits the discrete choices of a well: the side of the lower
bounding maximum (`eta_lo`), the lowest of several minima, and the first
crossing outward from the minimum. When two of these exchange order, the
derivative jumps (a kink) or, for the turning point of a level that crosses
the top of a secondary maximum inside the well, the value jumps. The lower
maximum switches sides on lines with almost equal maxima (LP QA, `mboz = 8`:
two lines with `(B_max,1 - B_max,2)/B = 4e-7`), and finite differences of
`Gamma_well` in `iota` then show this kink at steps above `1e-4`. As for the
other outputs, the Jacobian is the one-sided derivative of the frozen branch.

## Verification

- `dJ/dalpha` against a central difference of `J` over the label at fixed
  `eta`, with turning points found again by bisection (pyrabe test), and
  against `scipy` quadrature plus a label difference: `1e-7` relative.
- Fortran `TestSpectrumDerivative` (synthetic QH, QA), `TestSpectrumDerivativeQH/QA`
  (LP spectra) and `TestWoutChainQH/QA` (through libneo to the wout) check
  central differences per parameter block, adjoint duality and the wout chain
  for all seven outputs; best relative errors `1e-11` to `1e-7`.
- pyrabe (JAX) against `rabe_deriv.x` on the same spectra: values within
  `1e-13`, Jacobians within `1e-12` of their largest entry.
- `var_dJ` and `Gamma_s` against an independent scipy oracle (`quad` with
  the algebraic end weights, `brentq` turning points, 2D maximum of `B`) on
  the helical family with a synthetic radial block (`dB/ds = -0.02 - 0.06
  cos(theta + 4 phi) + 2 eps cos(theta)`, `iota' = 0.3`, `G' = 0.1`, `I =
  0.2`, `I' = 0.05`): `Gamma_s` within `3e-10`, `var_dJ` within `1.5e-7`
  for `eps = 1e-3 .. 3e-2` (pyrabe `test_alpha_criteria.py`). The oracle also
  checks `j_s` against a central difference of `J` in `s` with the shear term
  (`1.5e-9`), as does the pyrabe test (`1e-6`).
- Central differences (Fortran `TestSpectrumDerivative`, every parameter
  block including the radial one): `Gamma_s` `1e-11` to `3e-9`, `var_dJ`
  `2e-11` to `1.5e-8`; through the libneo wout chain with `ds_radial = 0.01`
  (`TestWoutChainQH/QA`) `3e-7` and `3e-9`. Fortran against pyrabe on the
  six stored spectra with a synthetic radial block: values within `5e-12`,
  Jacobians within `1.3e-11` of their largest entry.
