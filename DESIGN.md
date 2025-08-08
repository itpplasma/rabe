# CI/CD Pipeline Architecture Design

## Overview

This document outlines the architecture and implementation strategy for a robust CI/CD pipeline for the RABE project. The pipeline enforces strict Test-Driven Development (TDD) practices and comprehensive testing without shortcuts, stubs, or placeholders.

## Core Principles

1. **NO CHEATING**: No shortcuts, simplifications, stubs, or placeholders allowed
2. **TRUE TDD**: Always follow RED-GREEN-REFACTOR cycle with meaningful tests
3. **NON-TAUTOLOGICAL TESTS**: Tests must actually verify functionality, not just print "PASS"
4. **COMPREHENSIVE COVERAGE**: All tests labeled 'quick' and 'slow' must pass in CI
5. **QUALITY GATES**: Builds must pass all tests before merging

## Current Test Analysis

### Test Categories
- **quick**: 15 tests (~24 seconds) - Unit and integration tests
- **slow**: 3 tests (~64 seconds) - Performance-intensive tests  
- **plot**: Visualization tests (excluded from CI)
- **per_hand**: Manual tests (excluded from CI)

### Current Failures (5 tests failing)
1. TestGlobalMaximum (quick) - Numerical convergence issue
2. TestAntiSigmaAnalytic (quick) - Analytic comparison mismatch
3. TestPertAntiSigmaAnalytic (quick) - Perturbation calculation error
4. TestIIntegralHelicalAntiSigma (quick) - Integration accuracy issue
5. TestSetFieldlineLabels (quick) - Maximum resolution error

## CI Pipeline Architecture

### 1. GitHub Actions Workflow Structure

```yaml
name: CI
on:
  push:
    branches: [main, develop, feature/*]
  pull_request:
    branches: [main, develop]
```

### 2. Job Matrix Strategy

#### Build Matrix
- **Operating Systems**: ubuntu-latest, macos-latest
- **Build Types**: Debug, Release
- **Compiler**: gfortran
- **Python**: 3.11

#### Test Strategy
- **Standard Tests**: All quick tests on all OS/build combinations
- **Slow Tests**: Run on ubuntu-latest Release build only
- **Python Tests**: Run pytest suite if present

### 3. Dependencies

#### System Dependencies (Ubuntu)
- gfortran
- cmake (>= 3.16)
- libsuitesparse-dev
- liblapack-dev
- libblas-dev

#### System Dependencies (macOS)
- gcc (via Homebrew)
- cmake
- suite-sparse
- openblas
- lapack

#### Python Dependencies
- numpy
- matplotlib
- scipy
- pytest

### 4. Quality Gates

#### Code Quality Checks
1. **Line Length**: Strict 88 character limit enforcement
2. **No Commented Code**: Detection of commented-out code blocks
3. **No Dead Code**: Identification of unused functions/modules
4. **No Duplicate Code**: Detection of copy-paste violations

#### Test Quality Checks
1. **Non-Tautological**: Tests must have actual assertions
2. **Failure Detection**: Tests must fail when code is broken
3. **Coverage**: All public interfaces must be tested
4. **No Stubs**: No placeholder implementations allowed

### 5. Workflow Jobs

#### Job 1: Test Matrix
- Checkout code
- Setup environment
- Install dependencies
- Configure CMake
- Build project
- Run quick tests
- Run Python tests

#### Job 2: Slow Tests
- Separate job for performance tests
- 30-minute timeout
- Release build only
- Ubuntu only

#### Job 3: Code Quality
- Line length verification
- Commented code detection
- Code style compliance
- Test quality audit

## Implementation Backlog

### Epic 1: CI Pipeline Setup
- **Issue #1**: Create GitHub Actions workflow file
- **Issue #2**: Configure build matrix for OS/compiler combinations
- **Issue #3**: Set up dependency installation scripts
- **Issue #4**: Configure CMake build steps
- **Issue #5**: Integrate test execution with proper labels

### Epic 2: Fix Failing Tests
- **Issue #6**: Fix TestGlobalMaximum convergence
- **Issue #7**: Fix TestAntiSigmaAnalytic calculations
- **Issue #8**: Fix TestPertAntiSigmaAnalytic perturbation
- **Issue #9**: Fix TestIIntegralHelicalAntiSigma integration
- **Issue #10**: Fix TestSetFieldlineLabels resolution

### Epic 3: Test Quality Improvements
- **Issue #11**: Audit all tests for tautological assertions
- **Issue #12**: Remove any test stubs or placeholders
- **Issue #13**: Add missing edge case tests
- **Issue #14**: Improve test failure messages

### Epic 4: Code Quality Enforcement
- **Issue #15**: Add automated line length checks
- **Issue #16**: Add dead code detection
- **Issue #17**: Add duplicate code detection
- **Issue #18**: Add test coverage reporting

## TDD Implementation Process

For each issue:

### RED Phase
1. Write a failing test that verifies the requirement
2. Test must fail for the right reason
3. Test must be non-tautological (real assertions)
4. Commit the failing test

### GREEN Phase
1. Implement minimal code to pass the test
2. No shortcuts or stubs allowed
3. Implementation must be complete
4. All tests must pass
5. Commit the implementation

### REFACTOR Phase
1. Clean up implementation
2. Remove duplication
3. Improve naming and structure
4. Ensure all tests still pass
5. Commit the refactoring

## Success Criteria

1. All quick and slow tests pass consistently
2. CI pipeline runs on every push and PR
3. No tautological or stub tests remain
4. Code quality checks pass
5. Test coverage is meaningful and comprehensive
6. TDD process is followed strictly

## Anti-Patterns to Avoid

1. **Test Stubs**: Tests that always pass without checking anything
2. **Shallow Tests**: Tests that only check surface behavior
3. **Brittle Tests**: Tests that fail for unrelated changes
4. **Slow Feedback**: Tests that take too long to run
5. **Missing Tests**: Untested code paths or edge cases
6. **Commented Tests**: Tests that are disabled via comments

## Monitoring and Reporting

1. **Build Status Badge**: Display on README
2. **Test Results**: Published as artifacts
3. **Coverage Reports**: Track test coverage trends
4. **Performance Metrics**: Monitor test execution times
5. **Failure Analysis**: Track and analyze test failures

## Next Steps

1. Create detailed GitHub issues from this backlog
2. Prioritize fixing failing tests
3. Implement CI workflow incrementally
4. Enforce quality gates before merging
5. Document CI usage in CLAUDE.md