# trawl-fortran

Modern free-format Fortran translation of the computational code in the R package
`trawl` 0.2.2 by Almut E. D. Veraart.

The original package is GPL-3 and the translation is distributed under GPL-3.
The complete upstream source is retained in `upstream/trawl-0.2.2`.

## Implemented computational functionality

- Exponential, double-exponential, supIG, and long-memory trawl functions.
- Corresponding theoretical autocorrelation functions.
- Poisson and negative-binomial marginal moment estimators.
- Exponential, supIG, long-memory, and double-exponential trawl fitting from
  empirical autocorrelations.
- Generic and specialized trawl-set intersection calculations.
- Logarithmic-series moments, bivariate covariance/correlation, modified
  logarithmic-series moments, and bivariate modified moments.
- Bivariate negative-binomial simulation.
- Bivariate and trivariate logarithmic-series simulation.
- Univariate Poisson/negative-binomial trawl-process simulation.
- Fully-dependent and common-plus-idiosyncratic bivariate trawl simulation.

Plotting functions and plotting branches are intentionally omitted.

## Dependencies

The FPM package remains self-contained. The upstream R package imports
`DEoptim`, `rootSolve`, `Runuran`, `TSA`, `MASS`, `squash`, `ggplot2`, and
`ggpubr`. Only the first four participate in non-plotting computations.

- `DEoptim::DEoptim`: the supplied `DEoptim-fortran` 0.1.0 translation of
  DEoptim 2.2-8 is embedded in `src/` and used by the three GMM trawl fits.
  The complete dependency translation is retained under `vendor/`.
- `rootSolve::uniroot.all`: scalar root scanning plus bisection is implemented
  locally because this is the only rootSolve operation required by trawl.
- `Runuran::urlogarithmic`: logarithmic-series RNG is implemented locally.
- `TSA::acf`: the univariate demeaned sample ACF with lag zero omitted is
  implemented locally.

The DEoptim wrapper uses the same strategy/population/crossover/weight defaults
as upstream trawl. `set_trawl_seed()` also makes optimizer-backed fits
reproducible by deriving the standalone DEoptim seed from the trawl RNG.

## Build

```sh
fpm test
fpm run --example demo_trawl
```

The code uses `dp = kind(1.0d0)`, `implicit none`, free-form source, and no
nonstandard Fortran extensions in the public translation.
