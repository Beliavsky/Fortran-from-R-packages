# epo-fortran

A dependency-free modern Fortran implementation of the numerical algorithms in
`epo` 0.1.0.9000, the R package for Enhanced Portfolio Optimization (EPO).

The original package was written by Bernardo Reckziegel and distributed under
the MIT License. This project preserves that license and attribution.

## Features

- Simple Enhanced Portfolio Optimization.
- Anchored Enhanced Portfolio Optimization.
- Endogenous or exogenous risk-aversion scaling for anchored EPO.
- Correlation shrinkage from the sample correlation matrix toward identity.
- Optional full-investment normalization.
- Direct interfaces for return observations or a supplied covariance matrix.
- Sample covariance and covariance-to-correlation utilities.
- Typed results containing weights, covariance diagnostics, shrinkage matrices,
  the endogenous scaling coefficient, and error information.
- No external numerical-library dependency.

## Build with FPM

```text
fpm build
fpm test
fpm run epo_demo
fpm run --example simple_epo
fpm run --example anchored_epo
```

The return matrix is arranged as `(observations, assets)`.

## Basic use

```fortran
use epo, only : dp, epo_optimize, epo_result

real(dp) :: returns(100,4)
real(dp) :: signal(4)
real(dp) :: anchor(4)
type(epo_result) :: fit

! Fill returns, signal, and anchor.

fit = epo_optimize(returns, signal, 10.0_dp, 'anchored', 0.5_dp, &
  anchor=anchor, normalize=.true., endogenous=.true.)

if (.not. fit%ok) error stop trim(fit%message)
print *, fit%weights
```

A covariance matrix can be supplied directly:

```fortran
fit = epo_from_covariance(covariance, signal, 10.0_dp, 'simple', 0.5_dp)
```

## Mathematical conventions

For sample covariance matrix `C`, diagonal variance matrix `V`, and shrinkage
level `w`, the translated implementation forms

```text
C_tilde = (1 - w) C + w V.
```

This is algebraically identical to shrinking the correlation matrix toward the
identity and restoring the original standard deviations, as in the R source.

Simple EPO uses

```text
weights = (1 / lambda) inverse(C_tilde) signal.
```

Anchored EPO uses

```text
weights = inverse(C_tilde) [
  (1 - w) scale signal + w V anchor
].
```

`scale` is either `1 / lambda` or the endogenous coefficient described in the
upstream implementation. When normalization is enabled, the final vector is
divided by its sum. EPO does not impose long-only or box constraints, so
negative weights are possible.

## Input validation

The Fortran port requires:

- At least two observations when returns are supplied.
- Finite input values.
- Strictly positive asset variances.
- `0 <= w <= 1`.
- Positive `lambda` whenever it is used.
- A positive-definite shrunk covariance matrix.
- A nonzero weight sum when normalization is requested.

## License and provenance

- `LICENSE` contains the MIT License.
- `NOTICE` describes the derivative work.
- `original/epo-0.1.0.9000` contains the unmodified supplied source tree.
- `provenance/epo-main.zip` contains the supplied archive.
- SHA-256 manifests record original and translated files.

See `COVERAGE.md`, `PORTING_NOTES.md`, and `VALIDATION.md` for details.
