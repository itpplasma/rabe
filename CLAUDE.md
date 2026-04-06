# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

RABE (Redl's Analytical Bootstrap Estimation) is a scientific computing project for analytical bootstrap current estimation in stellarators, implemented primarily in modern Fortran with Python utilities for analysis and visualization.

## Build and Test Commands

### Building
- `make` - Build the project (Debug mode by default)
- `make CONFIG=Release` - Build in Release mode
- `make clean` - Clean build directory
- `make install` - Build and install

### Testing
- `make test` - Run standard tests (excludes slow/plot/current tests)
- `make test_slow` - Run performance-intensive tests
- `make test_all` - Run all tests except manual ones
- `make plot` - Run plotting/visualization tests
- `make current` - Run tests marked as current development focus

To run a single test: `ctest -R test_name --test-dir build`

### Python Testing
Python tests are located in `python/test/` and can be run with standard Python test runners.

## Architecture Overview

The codebase follows a layered architecture:

1. **Base Layer** (`src/utils/`)
   - `constants.f90` - Precision definitions and mathematical constants
   - `utils.f90` - Common utility functions

2. **Field Abstraction** (`src/`)
   - `field_base.f90` - Abstract interface for magnetic field calculations
   - `neo_field.f90` - Concrete implementation using NEO-2 splines

3. **Core Physics** (`src/`)
   - `fieldline_mod.f90` - Field line representation and properties
   - `make_fieldline.f90` - Field line construction and analysis
   - `fieldline_integrals.f90` - Integration along field lines
   - `deviation.f90` - Bootstrap current deviation factor calculations

4. **NEO Integration** (`src/neo/`)
   - Interfaces with NEO-2 code for magnetic field data
   - Handles `.bc` file reading and spline interpolation

The main application (`app/main.f90`) orchestrates these components to:
1. Read magnetic field configuration from `.bc` files
2. Create and analyze multiple field lines
3. Compute bootstrap current deviation factors
4. Output results to `rabe.out`

## Key Development Patterns

- **Derived Types**: Use `typename_t` naming convention (e.g., `fieldline_t`, `neo_field_t`)
- **Abstract Interfaces**: Field calculations use polymorphic `field_t` base type
- **Module Organization**: Each module has clear, focused responsibility
- **Test Categories**: Tests are labeled (slow, plot, current, per_hand) for selective execution

## Input/Output

- **Input**: Namelist file `rabe.in` containing configuration parameters
- **Field Data**: `.bc` files containing magnetic field spline data
- **Output**: Results written to `rabe.out`

## Plot Tests

See `test/plot/README.md` for the CMake pattern, available input files, myplot API reference, and examples.

## Important Notes

- Fortran arrays cannot be polymorphic - use wrapper types with abstract members
- Follow 88-character line limit and 4-space indentation
- The codebase integrates with NEO-2 modules for magnetic field calculations
- Build system enforces out-of-source builds (use `build/` directory)
