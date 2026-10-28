<img width="0" height="20" align="right">
<img src="assets/rabe_logo.png#gh-light-mode-only" alt="rabe logo by Georg Grassler" width="200" align="right"><img src="assets/rabe_logo_dark.png#gh-dark-mode-only" alt="rabe logo by Georg Grassler" width="200" align="right">

# rabe

[![Documentation](https://img.shields.io/badge/docs-FORD-blue)](https://itpplasma.github.io/rabe/)

An implementation of the nea**r**-omnigenous, **a**symptotic **b**ootstrap **e**xpressions of Ref. [1,2].

## Summary

In case of the Lorentz collision model,
the mono-energetic bootstrap and Ware pinch coefficients, $\bar D_{31}$
and $\bar D_{13}$, can be expressed
via a dimensionless geometry parameter $\lambda_{bB}$,

$$
\bar D_{31}=\bar D_{13} = \frac{1}{3} v \rho_L B \lambda_{bB}, \tag{1}
$$

scaling out particle velocity $v$, $\rho_L$ and $B$. This parameter can
be split into a ''symmetric'' contribution and a ''non-symmetric''  off-set,

$$
\lambda_{bB} = \lambda_\mathrm{sym} + \lambda_\mathrm{off}. \tag{2}
$$

For near-omnigenous stellarators at low plasma collisionality
$\nu_\ast=\pi R \nu_{\perp}/v$, where $\nu_\perp$ and $R$ are deflection frequency and
device major radius, respectively, the off-set can be expressed as

$$
    \lambda_{\text{off}} = \frac{\Lambda_\mathrm{A}}{\sqrt{\nu_\ast}} + \frac{\Lambda_\mathrm{B}}{\nu_\ast}, \tag{3}
$$

where $\Lambda_\mathrm{A}$ and $\Lambda_\mathrm{B}$ are the geometrical factors
due to the variation of the trapped-passing boundary layer width and the misalignment
of local maxima, respectively (see Ref. [1,2]).

Given a VMEC equilibrium file, `rabe` outputs those geometric coefficients
needed to evaluate $\lambda_\mathrm{off}$ at any collisionality.
Additionally, it can compute also the coefficients for the
''symmetric'' part

$$
\lambda_\mathrm{sym} = \lambda_{bB}^\mathrm{SC} + \lambda_{bB}^\mathrm{HGM}, \qquad \lambda_{bB}^\mathrm{HGM} = \Lambda_\mathrm{S}\sqrt{\nu_\ast} \tag{4}
$$

where the Shaing-Callen asymptotic $\lambda_{bB}^\mathrm{SC}$ is
approximated by the perfect omnigenous asymptotic
$\lambda_{bB}^\mathrm{SC} \rightarrow \lambda_{bB}^\mathrm{LC}$
obtained by Landreman and Catto [3] and
$\lambda_{bB}^\mathrm{HGM}$ is the contribution due to the finite
boundary correction derived by Helander, Geiger and Maasberg [4].

## Quick start

Requirements: gfortran and CMake >= 3.24; libneo, quadpack and Fortio are
fetched during the build (see [docs/building.md](docs/building.md), also for
macOS). The following runs the QH stellarator equilibrium of Ref. [5], which
ships with the repository as a test input:

```bash
make CONFIG=Release            # builds build/rabe.x
mkdir run_example && cd run_example
ln -s ../test/golden/input/rabe.in .
ln -s ../test/integration/vmec/input/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc .
../build/rabe.x
```

`rabe` reads `rabe.in` and the `field_file` named there from the working
directory and writes the coefficients per flux surface to `rabe.nc` (NetCDF)
and `rabe.dat` (plain text). To plot them:

```bash
python ../test/golden/plot_golden.py rabe.nc     # needs matplotlib and xarray
octave ../test/golden/plot_golden.m rabe.nc      # needs the Octave netcdf package
```

Both produce `rabe_output.png` with $\Lambda_\mathrm{A}$, $\Lambda_\mathrm{B}$,
$\Lambda_\mathrm{S}$ and, if `should_calc_shaing_callen = .true.`, the
Shaing-Callen asymptotic against $s_\mathrm{tor}$.

The core computation is also available from Python (`pip install rabe`;
see [docs/python.md](docs/python.md)).

## Documentation

- [docs/input-output.md](docs/input-output.md): the `rabe.in` namelist, the
  choice of helicity and the output variables.
- [docs/python.md](docs/python.md): Python package and examples.
- [docs/building.md](docs/building.md): build options, tests and dependencies.
- API documentation (FORD): <https://itpplasma.github.io/rabe/>

## Citing

If you use this code in any of your studies, please cite Refs [1,2].

See also [`CITATION.cff`](CITATION.cff) for machine-readable citation metadata.

## License

MIT License, see [LICENSE](LICENSE).

## References

[1] Georg S. Grassler, Christopher G. Albert, Sergei V. Kasilov, Winfried Kernbichler, *Asymptotic modeling of the bootstrap and Ware pinch effect in near omnigenous stellarators*, Phys. Plasmas; 33 (7): 070703. [doi:10.1063/5.0332431](https://doi.org/10.1063/5.0332431) (2026)

[2] C.G Albert et al., *On the convergence of bootstrap current to the Shaing–Callen limit in stellarators*, Journal of Plasma Physics, 91(3), p. E77. [doi:10.1017/S0022377825000200](https://doi.org/10.1017/S0022377825000200) (2025)

[3] M. Landreman and P. J. Catto, *Omnigenity as generalized quasisymmetry*, Phys. Plasmas 19, 056103 [doi.org/10.1063/1.3693187](https://doi.org/10.1063/1.3693187) (2012)

[4] P. Helander, J. Geiger, and H. Maassberg, “On the bootstrap current in
stellarators and tokamaks”, Phys. Plasmas 18, 092505 [doi.org/10.1063/1.3633940](https://doi.org/10.1063/1.3633940) (2011)

[5] M. Landreman et al., *Optimization of quasi-symmetric stellarators with self-consistent bootstrap current and energetic particle confinement*, Phys. Plasmas 29, [doi:10.1063/5.0098166](https://doi.org/10.1063/5.0098166) (2022)

[6] John R. Cary & Svetlana G. Shasharina, *Omnigenity and quasihelicity in helical plasma confinement systems, Phys. Plasmas 4*, 3323–3333, [doi:10.1063/1.872473](https://doi.org/10.1063/1.872473) (1997)
