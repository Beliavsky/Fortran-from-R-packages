# fgarch: experimental modern Fortran translation

This repository is an independent, experimental modern Fortran translation of the computational core of the R package **fGarch** version 4052.93.

## Status

**Experimental and incompletely validated.** The project builds, its examples run, and the included smoke/regression tests pass with gfortran 14. It has not been exhaustively compared with every R function, optimizer, distribution edge case, or historical fGarch result. Verify results independently before scientific, financial, or production use.

Plotting, graphical sliders, R S3/S4 presentation methods, R formula parsing, external data objects, GO-GARCH/ICA, and R-specific dependency plumbing are omitted.

## Implemented

- Variance-one normal, Student-t, GED, skew-normal, skew-Student-t, and skew-GED densities, CDFs, quantiles, and random generation.
- Corrected skew quantile threshold used by recent fGarch releases.
- Absolute moments for normal, Student-t, and GED distributions.
- GARCH/APARCH simulation with optional ARMA mean terms and arbitrary ARCH/GARCH orders.
- EGARCH simulation and filtering.
- Conditional log likelihood for the supported innovation distributions.
- Constrained GARCH(1,1) and APARCH(1,1) maximum-likelihood fitting.
- Standalone distribution fitting by maximum likelihood.
- Volatility forecasts, APARCH kappa, and persistence.
- Parametric VaR and expected shortfall.
- Jarque-Bera and Ljung-Box statistics.

## Build and run

```text
fpm build
fpm run
fpm test
fpm run --example fit_csv -- data/sample_returns.csv
```

Plain `fpm run` runs only `demo_fgarch`. The CSV example expects one numeric observation per line; a nonnumeric header is ignored.

Without fpm:

```text
make check
```

## Minimal use

```fortran
program example
   use fgarch, only : dp, garch_spec, garch_fit_result, make_garch_spec, &
      simulate_garch, fit_garch11, seed_rng, dist_std
   implicit none

   integer, parameter :: n = 1000
   type(garch_spec) :: spec
   type(garch_fit_result) :: fit
   real(dp) :: y(n), sigma(n), residuals(n)

   call seed_rng(1234)
   spec = make_garch_spec(1,1,cond_dist=dist_std)
   spec%omega = 1.0e-6_dp
   spec%alpha = [0.08_dp]
   spec%beta = [0.90_dp]
   spec%shape = 7.0_dp

   call simulate_garch(spec,n,y,sigma,residuals)
   fit = fit_garch11(y,cond_dist=dist_std,fit_shape=.true.)
   print *, fit%spec%alpha, fit%spec%beta, fit%spec%shape
end program example
```

## Important differences from R fGarch

- The public API is typed and array-based rather than S4/formula-based.
- The fitters currently estimate only GARCH(1,1) and APARCH(1,1); the lower-level filter and simulator accept arbitrary orders.
- Mean ARMA coefficients can be supplied for simulation/filtering but are not estimated by the current fitters.
- Optimization uses an internal Nelder-Mead implementation instead of R's optimizer collection.
- Standard errors and Hessian-based inference are not yet implemented.
- SNIG innovations are not implemented.
- VaR and expected shortfall are returned as lower-tail return quantiles, not automatically negated loss numbers.

See `API_MAP.md` and `VALIDATION.md` for exact coverage.

## License and provenance

The original package declares `GPL (>= 2)`. This translation is distributed under GPL-2.0-or-later and preserves attribution. See `LICENSE`, `NOTICE`, and `ORIGIN.md`.
