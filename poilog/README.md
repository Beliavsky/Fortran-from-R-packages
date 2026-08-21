# poilog-fortran

Modern Fortran/FPM translation of the computational core of the R package
`poilog` 0.4.2.1 (GPL-3), by Vidar Grotan and Steinar Engen.

## Implemented

- `dpoilog`: univariate Poisson-lognormal PMF.
- `dbipoilog`: bivariate Poisson-lognormal PMF.
- `rpoilog` and `rbipoilog`: random generation, including zero filtering and
  sampling conditional on the requested number of observed species.
- `poilog_mle_fit` and `bipoilog_mle_fit`: zero-truncated or untruncated maximum
  likelihood estimation.
- Parametric bootstrap support in both MLE routines.
- BFGS and Nelder-Mead optimization implemented in pure Fortran.

Plotting examples and R-specific data-frame/file-output machinery are omitted.

## Numerical translation notes

The original C implementation evaluates the PMFs by adaptive quadrature in
log-intensity space. This port retains the same high-mass interval search and
uses an embedded Gauss-Kronrod 15-point adaptive rule. For small univariate
counts, the original R wrapper also evaluates the standard-normal integral and
takes the larger result; this behavior is retained.

The bivariate implementation follows the original conditional-factorization
algorithm. The exact `rho = +/-1` boundary is handled explicitly as a
one-dimensional integral, avoiding the zero conditional variance that is
problematic in the upstream C routine.

Two obvious upstream R defects are corrected according to their intended
behavior:

1. `bipoilogMLE` checks `length(n1) == length(n2)` rather than comparing
   `length(n2)` to itself.
2. Bootstrap results are accepted on successful optimization rather than in the
   `try-error` branch.

## Build

```text
fpm build
fpm test
fpm run --example basic
```

The code uses only standard Fortran and has no external dependencies.

## API sketch

```fortran
use poilog, only : dp, dpoilog, dbipoilog, rpoilog, poilog_mle_fit, poilog_fit

p = dpoilog(3, 0.0_dp, 1.0_dp)
p2 = dbipoilog(2, 3, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.5_dp)
x = rpoilog(100, 0.0_dp, 1.0_dp, cond_s=.true.)
fit = poilog_mle_fit(x)
```

## License and attribution

The upstream package declares GPL-3. This translation is distributed under
GPL-3.0-only as a derivative work. The full GPLv3 text is in `LICENSE`.
Original package metadata and the source files used for the translation are
retained under `upstream_metadata/` for attribution and auditability.
