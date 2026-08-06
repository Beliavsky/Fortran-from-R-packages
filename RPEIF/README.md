# RPEIF modern Fortran

This package translates the computational core of the R package `RPEIF` to
modern Fortran and provides an FPM project.

The library computes influence functions for risk and performance measures:

- mean and robust mean
- standard deviation and semi-standard deviation
- value at risk and expected shortfall
- Sharpe, Sortino, and downside Sharpe ratios
- expected-shortfall and value-at-risk ratios
- Rachev ratio
- lower partial moments and Omega ratio

Plotting, `xts`/`zoo` indexing, R object dispatch, and other R-only
infrastructure are intentionally omitted.

## Build with FPM

```text
fpm test
fpm run --example influence_demo
```

Checked GNU Fortran flags can be used with:

```text
fpm test --flag "-std=f2018 -O0 -g -Wall -Wextra -Wpedantic -Wconversion-extra -fcheck=all -fbacktrace -ffree-line-length-none"
```

On systems without FPM, the included Makefile provides the same test suite:

```text
make MODE=checked test
make MODE=optimized test
```

The shell scripts in `scripts/` use FPM when available and otherwise use the
Makefile.

## Basic use

```fortran
use rpeif, only : dp, rpeif_options, influence_result, influence_series

real(dp) :: returns(100)
type(rpeif_options) :: options
type(influence_result) :: result

options%alpha = 0.05_dp
options%source_compatibility = .false.
call influence_series('ES', returns, result, options)
```

`result%x` contains the returns actually used after optional robust cleaning,
while `result%values` contains the influence-function series.

For evaluation at arbitrary points, use `influence_from_data`. For theoretical
shape evaluation under nuisance parameters, use `influence_from_nuisance` or
`evaluate_shape`.

## Estimator names

The dispatcher accepts:

```text
Mean SD SemiSD VaR ES SR SoR DSR ESratio VaRratio
RachevRatio robMean LPM OmegaRatio
```

Names are case-insensitive. `Omega`, `Rachev`, `RachR`, and several underscore
forms are accepted as aliases.

## Options

The `rpeif_options` type controls tail probabilities, risk-free rate, lower
partial-moment threshold and order, Sortino threshold, robust family and
efficiency, robust cleaning, AR prewhitening, and source compatibility.

The default `source_compatibility=.true.` preserves the formulas and empirical
code paths in the upstream R package. Setting it to `.false.` corrects three
apparent upstream inconsistencies:

1. empirical UPM sign in the Omega influence function;
2. use of the quantile rather than density at the quantile in empirical
   VaR-ratio influence values;
3. omission of the risk-free rate from the scale term of the Sharpe-ratio
   influence function.

See `PORTING_NOTES.md` for details.

## Dependency

RPEIF uses a small vendored subset of the completed RobStatTM Fortran
translation for rho/psi functions and efficiency tuning. The RPEIF-specific
robust location and cleaning iteration is implemented in `rpeif_robust.f90`.
No BLAS or LAPACK library is required.

## License

The upstream RPEIF package is GPL (>= 2). The reused RobStatTM translation is
GPL-3.0-or-later. The combined Fortran package is therefore distributed under
GPL-3.0-or-later. Complete license texts and provenance are included.
