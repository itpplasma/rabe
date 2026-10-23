# Boundary gradients of rabe's offset (VMEC++ -> Boozer -> rabe)

`rabe_boundary.BoundaryObjective` gives `objective(boundary)` and the exact
`grad(boundary)` of rabe's offset coefficients with respect to the VMEC++
fixed boundary `boundary = [rbc, zbs]` (shape `(2, mpol, 2 ntor + 1)`, VMEC++
layout). The default objective is the 1/nu offset summed over the chosen
surfaces,

    J = sum_s w_s L_s,   L = Lambda_A / sqrt(nu_star) + Lambda_B / nu_star,

as `pyrabe.simsopt.RabeOptimizable.Lambda_offset`, optionally squared
(`squared=True`). Any smooth function of the `(S, 5)` outputs `(Lambda_A,
Lambda_B, Lambda_S, nu_star_crit, eta_b)` can be passed as `reduce`, and
`vjp(boundary, cotangent)` pulls back an arbitrary output cotangent.

## Chain

1. VMEC++ (PyPI `vmecpp>=0.9.1` with the Enzyme force Jacobian,
   `vmecpp.has_exact_force_jacobian()`) solves the equilibrium with the m=1
   gauge pinned. The gradient uses its implicit boundary adjoint (reverse mode).
2. Only the primary wout arrays enter: `rmnc`, `zmns`, half-mesh `lmns`,
   `iotaf`, `phi[-1]` and `Rmajor_p`. The derived `bmnc`, `bsub*`, `gmnc` do
   not.
3. [pyrabe](https://github.com/itpplasma/pyrabe) `WoutChain`: the libneo Boozer
   transform with the `-I*lambda` term of libneo PR 450 in JAX (including
   `|grad s|`, `psi_edge` and `R`), then rabe's derivative mode (PR 125) in JAX.
   pyrabe's derivative mode agrees with `rabe_deriv.x` to 1e-10.

One gradient costs one VMEC++ solve, one pyrabe VJP per surface and one
adjoint solve. `objective` and `grad` at the same boundary share the solve.
rabe's discrete choices (`p/q`, labels, bounding maxima) are recomputed at
every new boundary and frozen within it, so the gradient is that of the smooth
branch through the current point.

Environment: `uv pip install "vmecpp==0.9.1" -e <pyrabe checkout>` (pyrabe main
a1f7bfc or later), then run the scripts from this directory.

## Scripts

- `check_gradient.py`: full-chain central differences (every point re-solves
  VMEC++) against the gradient, per output and surface, frozen and recomputed
  structure.
- `ns_scan.py`: values, gradients and cost per gradient against `ns`.
- `smoke_descent.py`: a few steepest-descent steps with backtracking.
- `test_rabe_boundary.py`: pytest (`RABE_OPT_INPUT=<input> pytest opt`).

The force tolerance matters: at `ftol = 1e-14` the half-mesh `lmns` is
converged only to about 1e-3, and the gradient moved by 3 % between
`ftol = 1e-14` and `1e-17` (LP QH, ns = 25). The scripts use `ftol = 1e-20`.

## Limits

Stellarator-symmetric fixed-boundary inputs only. `phiedge` and the profiles
are held fixed (VMEC++ also differentiates power-series profiles, not wired
here). The iota direction of rabe uses the frozen rational surrogate `p/q`.
Derivatives are of the discretization (VMEC++ at the given `ns`, the libneo
splines, rabe's Gauss-Legendre rule), not of the continuum.
