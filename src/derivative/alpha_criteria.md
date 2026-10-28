# Alpha-particle criteria of the derivative mode: `Lambda_alpha` and `Gamma_well`

The derivative mode returns two outputs in addition to the offset coefficients,
`Lambda_alpha` (index 6) and `Gamma_well` (index 7), with exact gradients like
the others. They implement two of the smooth optimization criteria that S. V.
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
