# Building and testing rabe

## Prerequisites

- Fortran compiler (gfortran)
- CMake >= 3.24

On Debian/Ubuntu:

```bash
sudo apt-get install gfortran cmake
```

On macOS (Homebrew):

```bash
brew install gfortran
export FC=$(ls $(brew --prefix)/bin/gfortran-* | sort -V | tail -1)
```

Make sure that `FC` points to your installed gfortran compiler. If
[Ninja](https://ninja-build.org) is available, we recommend it as the CMake
generator:

```bash
export CMAKE_GENERATOR=Ninja
```

## Build

```bash
make                  # Debug build (default CONFIG=Debug)
make CONFIG=Release   # optimized build
make clean
```

The executable `rabe.x` is written to `build/`. `make python` configures a
Release build with the Python bindings (`-DBUILD_PYTHON_BINDINGS=ON`); for
installing the Python package see [python.md](python.md).

## Tests

```bash
make test         # main test suite (CTest label "quick")
make test_all     # also the long-running tests (labels "quick" and "slow")
make test_slow    # only the long-running tests
make test_failed  # rerun failed tests verbosely
make golden       # golden-record comparison
```

The tests and their descriptions are in `test/`.

## Overriding fetched dependencies

By default the build fetches pinned commits of libneo and Fortio (see
`cmake/libneo.cmake` and `cmake/fortio.cmake`). Other refs can be selected
explicitly, for example for a compatibility build. The `Makefile` forwards
these variables only when they are given on the command line.

- `-DFORTIO_REF=<branch|tag|sha>` selects a different Fortio ref:

  ```bash
  make FORTIO_REF=main
  # or directly: cmake -S . -B build -DFORTIO_REF=main
  ```

- `-DLIBNEO_REF=<branch|tag|sha>` selects a different libneo ref:

  ```bash
  make LIBNEO_REF=main
  # or directly: cmake -S . -B build -DLIBNEO_REF=main
  ```

- `-DLIBNEO_PATH=<dir>` uses a local libneo checkout instead of fetching:

  ```bash
  make LIBNEO_PATH=/path/to/libneo
  # or directly: cmake -S . -B build -DLIBNEO_PATH=/path/to/libneo
  ```

## Third-party code

Fetched automatically during the build:

- [libneo](https://github.com/itpplasma/libneo) for field file I/O (MIT)
- [quadpack](https://github.com/jacobwilliams/quadpack) for numerical
  integration (BSD-3-Clause)
- [fortio](https://github.com/lazy-fortran/fortio) for NetCDF input/output (MIT)

[pyplot-fortran](https://github.com/jacobwilliams/pyplot-fortran) is used
optionally for visualization; its source and license are in `test/plot_lib`
(BSD-3-Clause).
