The SD2.2 pointwise kernels are generated with FortSym. Ordinary RABE builds
compile the checked-in output and do not require a computer algebra system.
The generator verifies the first three series-recurrence coefficients of
`i*a*gamma=(b*gamma')'` and the fourth power of the positive `A_m` amplitude
against its primitive FEdef normalization, and differentiates the exact
logarithmic action to recover `sqrt(-log(x))`. These are exact symbolic identities,
not numerical probes. Numerical physics checks remain separate CTest tests.

With an existing FortSym build, run from the RABE repository root (replace
`FORTSYM_BUILD` with that build directory):

```sh
export FORTSYM_BUILD=/path/to/fortsym/build
gfortran -I "$FORTSYM_BUILD/include" gen/sqrt_nu_kernel.f90 \
  -L "$FORTSYM_BUILD/lib" -lfortsym \
  -Wl,-rpath,"$FORTSYM_BUILD/lib" -o /tmp/rabe-gen-sqrt-nu
/tmp/rabe-gen-sqrt-nu
git diff --exit-code -- src/generated/sqrt_nu_kernel.f90
```

The generator uses FortSym's SymEngine backend and fails on an unknown exact
identity. Its output implements the real decaying WKB profile and derivative,
the logarithmic FEdef quadrature integrand and closed-form diagnostic, and
`A_m` from `nearomnigen.tex` Eq. `largeparam`. The recurrence verification
checks the note's expression containing `b3` in the cubic coefficient.

FortSym currently fixes its emitter line limit at 100 columns and exposes no
width option in `kernel_spec_t`. A few generated declarations and expressions
therefore exceed RABE's preferred 88 columns, while remaining valid Fortran.
The generated output is kept intact.
