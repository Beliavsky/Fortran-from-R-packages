# tseries-modern-fortran

Experimental modern Fortran translation of the computational portions of the R package `tseries` 0.10-62.

## Status

**Experimental.** The project builds and its included smoke and extended tests pass with gfortran 14 under runtime checking. It has not been extensively validated against every result, edge case, optimizer behavior, or platform supported by the original R package. Independently verify numerical results before production, scientific, financial, medical, or safety-critical use.

The translation intentionally omits plotting, internet data download, R S3 presentation methods, and most irregular-time-series/date/file-format plumbing.

## Implemented computational API

- Conditional-sum-of-squares ARMA fitting
- Gaussian GARCH fitting and conditional-variance prediction
- Augmented Dickey-Fuller test
- Phillips-Perron test
- KPSS stationarity test
- Phillips-Ouliaris cointegration test
- BDS test using a direct low-memory algorithm
- Jarque-Bera test
- Runs test
- Terasvirta neural-network nonlinearity test
- White neural-network nonlinearity test
- Stationary and fixed-block bootstrap
- Random-permutation, Fourier, and amplitude-adjusted surrogates
- Quadratic/logistic map generation
- Maximum drawdown, Sharpe ratio, and Sterling ratio
- Mean-variance portfolio optimization

See `API_MAP.md` for function-by-function coverage.

## Build and run

```text
fpm build
fpm run
fpm test
```

Plain `fpm run` executes only `demo_tseries`; the CSV utility is an example target, so it does not fail from a missing command-line argument.

Run the CSV example with:

```text
fpm run --example fit_csv -- example/sample.csv
```

Depending on the installed fpm version, this equivalent form may also work:

```text
fpm run fit_csv -- example/sample.csv
```

The CSV example expects one numeric observation per row. Non-numeric header lines are skipped.

## Public module

```fortran
use tseries, only : dp, arma_result, arma_fit

real(dp) :: x(100)
type(arma_result) :: fit

! Fill x.
fit = arma_fit(x, p=1, q=1)
print *, fit%coefficients
```

All public procedures are exported by the `tseries` module. More specialized implementation modules remain available but should be regarded as internal.

## Compiler settings

The manifest enforces modern source rules:

```toml
[fortran]
implicit-typing = false
implicit-external = false
source-form = "free"
```

The included Makefile uses:

```text
-std=f2018 -Wall -Wextra -Wimplicit-interface -fcheck=all -fbacktrace
```

## Numerical differences from R

This is not a byte-for-byte language conversion. Important differences include:

- ARMA and GARCH use an internal Nelder-Mead optimizer.
- GARCH parameters are transformed to enforce positivity and persistence below one.
- Portfolio optimization uses direct equality-constrained systems and a simple active-set reduction for no-short portfolios, rather than `quadprog::solve.QP`.
- BDS uses a simpler direct pair-counting implementation instead of the original bitmap accelerator.
- Fourier transforms use a dependency-free O(n^2) DFT, suitable for examples and moderate samples but not large production workloads.
- White's test uses a deterministic seed by default so examples are reproducible.
- R object classes, formula parsing, names, time-series attributes, warnings, and print/plot methods are not reproduced.

## Repository contents

```text
src/        library modules
app/        no-argument demonstration used by fpm run
example/    command-line CSV example and sample data
test/       smoke and extended tests
reference/  provenance notes and room for future R reference outputs
```

## License and provenance

The original R package declares `GPL-2 | GPL-3`. This translation retains the original notices and is distributed under GPL version 2 only or GPL version 3 only, at the recipient's option. See `LICENSE`, `LICENSE-GPL-2`, `LICENSE-GPL-3`, `NOTICE`, and `ORIGIN.md`.

This is an independent and unofficial translation. It is not endorsed by the original package authors or CRAN.
