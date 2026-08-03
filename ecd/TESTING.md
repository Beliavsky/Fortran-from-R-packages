# Testing

The release is tested in two configurations.

## Strict configuration

```text
-std=f2018
-Wall -Wextra -Werror -Wno-compare-reals
-fcheck=all
-ffpe-trap=invalid,zero,overflow
-fbacktrace
-O0 -g
```

## Release configuration

```text
-std=f2018
-Wall -Wextra -Werror -Wno-compare-reals
-O3
```

## Permanent tests

- `test_ecd_core.f90`
  - standard cusp normalization and moments
  - PDF/CDF/quantile identities
  - ECLD constants and independent fixed references
- `test_extended_models.f90`
  - scaled quartic error functions
  - OGF-star analytic identities
  - incomplete moments
  - quartic MGF/OGF put-call parity
  - quartic Q and RN0 constants
- `test_processes.f90`
  - Laplace and Lihn-Laplace distributions
  - stable-count fixed references and inversion
  - SLD/QSLD characteristic functions and cumulants
  - deterministic random generation checks
- `test_options.f90`
  - Black-Scholes fixed references
  - implied-volatility inversion
  - split polynomial option fitting
- `test_utilities_lamp.f90`
  - sample and lag calculations
  - R type-7 empirical quantiles
  - tail partitioning
  - LAMP simulation identities
- `test_fitting.f90`
  - bounded Nelder-Mead
  - ECD constant estimation
  - fitting-result structure and finite likelihoods

All application and example targets are also compiled and executed.

## Independent references

Reference values were generated independently with Python/SciPy formulas or
closed-form identities. Important fixed references include:

- ECLD lambda=3, sigma=0.4 at x=0.2
- standardized Lihn-Laplace density
- stable-count gamma special case
- Black-Scholes call prices and implied-volatility inversion
- standard cusp normalization, variance, and kurtosis
- quartic OGF-star and put-call parity identities

See `REFERENCE_GENERATION.md` for details.

## GNU linker note

GNU Fortran can implement procedure callbacks with stack trampolines. On Linux,
GNU `ld` can consequently print an executable-stack warning for objects using
internal callback procedures. The strict and release executables run correctly;
this warning is a compiler implementation detail rather than a failed test.
