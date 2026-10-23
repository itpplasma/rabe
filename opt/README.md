# Boundary gradients of rabe's offset (VMEC++ -> Boozer -> rabe)

`rabe_boundary.BoundaryObjective` gives `objective(boundary)` and the exact
`grad(boundary)` of rabe's offset coefficients with respect to the VMEC++
fixed boundary `boundary = [rbc, zbs]` (shape `(2, mpol, 2 ntor + 1)`, VMEC++
layout). The default objective is the 1/nu offset summed over the chosen
surfaces,

    J = sum_s w_s L_s,   L = Lambda_A / sqrt(nu_star) + Lambda_B / nu_star,

as `pyrabe.simsopt.RabeOptimizable.Lambda_offset`, optionally squared
(`squared=True`). Any smooth function of the `(S, 7)` outputs `(Lambda_A,
Lambda_B, Lambda_S, nu_star_crit, eta_b, Lambda_alpha, Gamma_well)` can be
passed as `reduce`, and `vjp(boundary, cotangent)` pulls back an arbitrary
output cotangent.

## Alpha objective

`alpha_reduce(AlphaWeights(...))` is the composite objective of Kasilov's
criteria (see [alpha_criteria.md](../src/derivative/alpha_criteria.md)),

    F = sum_s w_s [ w1 nu*_crit + w2 (Lambda_A^2 + Lambda_B^2) + w2b Lambda_alpha^2 + w3 Gamma_well ],

with every term divided by a reference (`AlphaWeights().normalized(values)`
takes the values at the start boundary, so that the weights are relative).
The offsets enter squared because their plain sum is unbounded below.
`constraints=[...]` adds penalties `c(boundary, harm)`, JAX functions of the
boundary and of the primary wout arrays whose gradient joins the same VMEC++
adjoint: `aspect_penalty` and `volume_penalty` from the boundary geometry
(`boundary_geometry`, equal to VMEC++'s `aspect` and `volume_p` to 1e-15 on LP
QA and QH) and `iota_penalty` on `iotaf`. A QI mirror-ratio constraint needs
the Boozer spectrum per surface and is not wired yet.

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
- `alpha_optimize.py`: L-BFGS-B on the boundary for the experiments E1
  (`nu*_crit`, `Lambda_A^2 + Lambda_B^2`) and E2 (`nu*_crit`, `Lambda_alpha^2`,
  `Gamma_well`), with aspect, volume, iota-range and Boozer mirror-ratio
  penalties (`iota_range_penalty`, `MirrorRatioPenalty`) around the start
  point. It checkpoints every iteration, restarts from the checkpoint, stages
  and submits SIMPLE alpha-loss runs of the start, every N-th and the final
  iterate, and has an `ftol-check` mode that picks the VMEC++ force tolerance.
- `vmecpp_jacobian.py`: the traced `vmecpp.run` with hot restarts from the
  last converged solve and one factorization of the implicit adjoint per
  equilibrium, shared by all cotangents of a `jax.vmap`ped pullback (private
  vmecpp 0.9 internals, version checked on import).
- `lsq_opt.py`: E1, E2 and E3 (`nu*_crit`, `var_dJ`, `Gamma_s`, needs a pyrabe
  with those outputs) as a residual vector with the constraints as penalty or
  augmented-Lagrangian residuals, the exact Jacobian, Levenberg-Marquardt with
  a step cap and an augmented Lagrangian.
- `run_opt.py`: the checkpointed campaign entry point
  (`--engine vmecpp --objective E1|E2|E3 --solver lm|trf|lbfgs|al`) with
  continuation in the free modes and in `ns`, multi-start and SIMPLE staging.
- `check_jacobian.py`: FD V-curve of the residual Jacobian, with the rabe
  structure re-found and frozen; `test_lsq_opt.py`: the batched adjoint
  against vmecpp's own VJP and the Jacobian against a frozen-structure FD.
- `test_rabe_boundary.py`: pytest (`RABE_OPT_INPUT=<input> pytest opt`), the
  offset gradient and the alpha objective with aspect and iota penalties
  against re-solved central differences.

The force tolerance matters: at `ftol = 1e-14` the half-mesh `lmns` is
converged only to about 1e-3, and the gradient moved by 3 % between
`ftol = 1e-14` and `1e-17` (LP QH, ns = 25). The scripts use `ftol = 1e-20`.

## Limits

Stellarator-symmetric fixed-boundary inputs only. `phiedge` and the profiles
are held fixed (VMEC++ also differentiates power-series profiles, not wired
here). The iota direction of rabe uses the frozen rational surrogate: within
one structure the field-line iota is `p/q + (iota - iota_input)`
(`pyrabe.spectrum_derivative.locate`), whereas a fresh structure resets it to
`p/q`. The objective with the structure recomputed at every point therefore
differs from the differentiated branch by a smooth term, 6e-4 of the
directional derivative on LP QA (ns = 25) and 6e-7 on LP QH for the offsets;
for the least-squares residuals on LP QH (ns = 25) it is 8e-4 (E1) and 2e-3
(E2) of the largest entry of `J d`, against 3e-5 for the frozen structure.
Derivatives are of the discretization (VMEC++ at the given `ns`, the libneo
splines, rabe's Gauss-Legendre rule), not of the continuum.
