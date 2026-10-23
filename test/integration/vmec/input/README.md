# Test data

The test suite uses the following VMEC equilibrium files

- **Landreman-Paul QA** (`wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc`): A quasi-axisymmetric stellarator equilibrium at reactor scale from Ref. [1]. Taken from the test data of [simsopt](https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc).
- **Landreman-Paul QH** (`wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc`): A quasi-helically symmetric stellarator equilibrium at reactor scale from Ref. [1]. Taken from the test data of [simsopt](https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc).

## Refined QH fixture

`wout_LandremanPaul2021_QH_reactorScale_ns101_reference.nc` is a public-derived
radial refinement of the same reactor-scale QH equilibrium. Native tests copy
this tracked file locally and require neither a network connection nor VMEC++.
The original ns50 file above remains unchanged. The NEO reference values in
`test_against_neo_qh.f90` were regenerated for this ns101 fixture with the same
booz_xform, `.bc` writer and `neo_field_t` route that reproduces the original
ns50 values, as recorded in that test.

The input comes from [simsopt revision
81436b4f7c655a4099eadadb168b7cc27923e555](https://github.com/hiddenSymmetries/simsopt/blob/81436b4f7c655a4099eadadb168b7cc27923e555/tests/test_files/input.LandremanPaul2021_QH_reactorScale_lowres).
Its SHA256 is f1611722b94ac7d40e30c7f51be273d0db50de19bd3ba4cc04bbd512d3169e0d.
Only three solver arrays change in the adjacent public input:

| Setting | Original | Refined |
|---|---|---|
| NS_ARRAY | 12, 50 | 12, 101 |
| FTOL_ARRAY | 1e-17, 1e-15 | 1e-17, 1e-18 |
| NITER_ARRAY | 1200, 5000 | 1200, 20000 |

The refined input SHA256 is
1ccbdcdd29d79f13b02a49b3f139d628e556b036a41ea592236162f016327fde.
Its boundary, NFP=4, MPOL=8, NTOR=8, fixed-boundary stellarator symmetry,
zero pressure and prescribed-current profiles, and edge toroidal flux
41.86387773208 Wb retain the public input values. Compared with the original
wout, the largest outer-boundary coefficient differences are 2.22e-16 m in R
and 5.56e-17 m in Z. The residual total current is 7.93e-9 A.

The tracked wout contains 2,623,464 bytes and has SHA256
a8fb3ac3c0cc59630b53ebc495b4505e4754660f9d17e30bd82ad18bc2f11728.
The saved output has NS=101, IER_FLAG=0, FTOLV=1e-18 and force residuals
(fsqr, fsqz, fsql) = (9.826771929772847e-19, 9.845205945360565e-19,
1.1152947042899614e-19). Each residual satisfies the fixed finite,
nonnegative, <=1e-17 qualification bound. All unmasked numerical data are
finite; the only absent masked scalar is extcur, also absent in the original.

The unchanged libneo HYBRID quality check at d03789a1fef15e188752614328dd50ea9854919f
(physics identical at e2f39a76645bed0919bca9fc239a98114a930631) returned status 0
with allow_poor_equilibrium=0, 16 checks and rho in [0.2, 0.95]. Its six values
were (1.7584751822873391e-6, 3.198484870264345e-4, 5.326484766293867e-4,
1.3339712208247181e-4, 3.8124851715889043e-7, 0), below the default bounds
(1e-5, 1e-3, 0.05, 6e-4, 1e-5, 2e-4). The original ns50 wout fails the
force-balance bound under the new default construction.

### Optional regeneration

Generation used the installed VMEC++ 0.9.1 wheel/API with Python 3.12.14,
NumPy 2.5.4, SciPy 1.18.1 and netCDF4 1.7.5 on Linux x86-64-v3, with one
thread. VMEC++ release commit 69a51e06331a89d0b974cbcd6e80c3edb995af3e
identifies its release source, rather than attesting the installed binary build.
The actual runtime components had these SHA256 identities:

| Component | SHA256 |
|---|---|
| Python 3.12.14 executable | 8f9781a98200d9ecda7e00464e4c64b1327abae788ae8e6979d5c859311410c7 |
| VMEC++ Python binding | d0bf7651d1d38be3316728ffa3f7521c12bf6b3801caad642a7b1cb11d0a9cf3 |
| Selected x86-64-v3 core | 642ff045bea5da0bbc9d8e55ce7e73ee2def23f1957f0e1f9922931f9405f355 |
| indata2json converter | 3744aa47e915f869894f3f8ee66fdea9038ca562d9058ccf977cf1aaa9d415fe |
| libc.so.6 | 6b4a45352fd0c540a9c7c718f35ce8c8e46a4e482f9d3885a910c32d1a0e1421 |

The adjacent recipe uses one vmecpp.run call with the two specified radial
stages. Run it from this directory with a fresh output directory outside Git:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python regenerate_qh_ns101.py \
  input.LandremanPaul2021_QH_reactorScale_ns101 /path/to/new-output
```

The classic-input converter warns that PT_TYPE and PH_TYPE are unknown. It
serializes PMASS_TYPE='power_series', AM=[0], PCURR_TYPE='power_series', AC=[0],
GAMMA=0, NCURR=1 and CURTOR=0. PT/PH select anisotropy/hot-pressure or flow
profiles in other VMEC variants and do not enter this isotropic VMEC++ model,
as defined by the [converter serialization](https://github.com/jonathanschilling/indata2json/blob/a29c5e3e932df15dd2beab0341b1d7bbe7eb32eb/src/indata2json.f90#L267)
and [VMEC profile functions](https://github.com/PrincetonUniversity/STELLOPT/blob/8060f5e5b1bfe11b2f8809c8dfa90459ee72f9ba/LIBSTELL/Sources/Miscel/profile_functions.f#L874).
The saved pressure and mass profiles and input AM/AC/AI arrays are zero, so
these ignored default selectors do not change this equilibrium. This does not
establish equivalence for nonzero anisotropic or flow profiles.

Regenerated output must satisfy the saved-residual, physical-identity and
unchanged field-quality contracts above before replacing the tracked fixture.
Different runtime builds may produce different file bytes or floating-point
roundoff, and any replacement needs its own recorded checksum and review.

## References

- [1] M. Landreman & E.J. Paul, *Magnetic fields with precise quasisymmetry for plasma confinement*, Phys. Rev. Lett. 128, 035001, [doi:10.1103/PhysRevLett.128.035001](https://doi.org/10.1103/PhysRevLett.128.035001) (2022)
