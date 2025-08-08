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

## Input/Output

- **Input**: Namelist file `rabe.in` containing configuration parameters
- **Field Data**: `.bc` files containing magnetic field spline data
- **Output**: Results written to `rabe.out`

## Development Principles

- Strictly follow TDD, SOLID, KISS, SRP and DRY as strict requirements in the development process.
- !!! You MUST always write tests first. Strictly follow Test Driven Development as defined below !!!
- Always keep functions, subroutines and modules small and manageable.
- Always keep work units small and manageable.
- Always work sequentially. Never do multiple things in parallel.
- Never be lazy. Don't take shortcuts. Don't stub code. Don't write placeholders. Fulfill every task fully to 100%
- Delete things that are obsolete if they are under version control. Don't hamster old messy stuff.
- Always write code that speaks for itself and avoid comments where possible.
- Never leave commented out sections of code in source files.
- You have a strict limit of 88 characters per line. 
- You indent with 4 spaces.
- When you work, you carefully go through each task in order, check of completed items, and update the plan when needed. Always plan before you implement. Always write tests first.
- Usually, you don't need to clean the build output before build or test. Do this only in case of persisting build problems that you cannot solve otherwise.
- Use typename_t as convention for Fortran derived types
- Always prefer simple and elegant solutions over complex ones.
- Fortran knows no polymorphic arrays. You need a concrete wrapper type with one abstract type member variable and build an array of this wrapper type objects.
- Never nest deeper than 3 levels in loops, if blocks, etc.
- To avoid dummy argument warnings in Fortran, you can place them in an empty associate block
- Extend Fortran arrays with arr = [arr, new_element] syntax. Be sure that new_element is just a variable and not an expression (one may need a temporary)
- I hate having several versions of the same thing lying around like a Python programming rookie who has never used git. Clean up redundant code immediately.
- NEVER keep old versions and start new files with updates or rewrites. ALWAYS work on existing files in place.
- NEVER use phrases like "for now assume" or take methodological shortcuts.
- When cleaning the build, use fpm clean --skip to avoid removing files that should be kept.

## Test Driven Development

Follow the Red-Green-Refactor cycle:

1. **RED**: Write a test that will fail until we wrote our code
2. **GREEN**: Implement actual code and pass this test  
3. **REFACTOR**: Clean up code while keeping all tests green

## Rules
- Never write code without a failing test
- Write only enough code to pass the test
- Run tests after each step
- A REAL test must:
  - Actually exercise the functionality being tested
  - Fail when the functionality is broken
  - Pass only when the functionality works correctly
  - Not just print "PASS" or "FAIL" messages
- For a compiler, real tests compile actual Fortran code and verify the output

## Personal Information

- My github username is GeorgGrassler

## Communication Guidelines

- Be concise and professional. No celebrations, no exaggerations.
- Don't declare victory when fixing minor infrastructure issues if the main functionality still doesn't work
- A test is NOT just code that prints "PASS" - it must actually verify functionality
- Understand the difference between infrastructure fixes and actual feature implementation

## Important Notes

- Fortran arrays cannot be polymorphic - use wrapper types with abstract members
- Follow 88-character line limit and 4-space indentation
- The codebase integrates with NEO-2 modules for magnetic field calculations
- Build system enforces out-of-source builds (use `build/` directory)
