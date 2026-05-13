# rabe

An implementation of the nea**r**-omnigenous, **a**symptotic **b**ootstrap **e**xpressions of Ref. [1].

## Build

use

```bash
make
```

to build the executable and

```bash
make test
```

or

```bash
make test_all
```

to execute tests.

## Third Party

- [`pyplot-fortran`](https://github.com/jacobwilliams/pyplot-fortran) for visualization (source and license under `plot`, BSD-3-Clause)

- [`quadpack`](https://github.com/jacobwilliams/quadpack) for numerical integration (fetched during build, BSD-3-Clause)

## References

[1] C.G Albert et al., *On the convergence of bootstrap current to the Shaing–Callen limit in stellarators*, Journal of Plasma Physics, 91(3), p. E77. [doi:10.1017/S0022377825000200](https://doi.org/10.1017/S0022377825000200) (2025)

[2] M. Landreman et al., *Optimization of quasi-symmetric stellarators with self-consistent bootstrap current and energetic particle confinement*, Phys. Plasmas 29, [doi:10.1063/5.0098166](https://doi.org/10.1063/5.0098166) (2022)

[3] A. Redl et al., *A new set of analytical formulae for the computation of the bootstrap current and the neoclassical conductivity in tokamaks*, Phys. Plasmas 28, [doi:10.1063/5.0012664](https://doi.org/10.1063/5.0012664) (2021)

[4] O. Sauter et al., *Neoclassical conductivity and bootstrap current formulas for general axisymmetric equilibria and arbitrary collisionality regime*, Phys. Plasma 6, [doi:10.1063/1.873240](https://doi.org/10.1063/1.873240) (1999)

[5] John R. Cary & Svetlana G. Shasharina, *Omnigenity and quasihelicity in helical plasma confinement systems, Phys. Plasmas 4*, 3323–3333, [doi:10.1063/1.872473](https://doi.org/10.1063/1.872473) (1997)
