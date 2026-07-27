# betategarch-modern-fortran

A modern Fortran translation of the computational core of the R package
`betategarch` 3.4 by Genaro Sucarrat. The attached R package describes itself
as providing simulation, estimation, and forecasting for first-order
Beta-Skew-t-EGARCH models.

## License

The original package declares `License: GPL-2`. This project therefore uses
**GPL-2.0-only**. The complete license is in `LICENSE`, and every Fortran source
file carries an SPDX identifier, the original copyright attribution, and an
explicit GPL version 2-only notice.

## Implemented and tested

- Fernandez-Steel skew-t random generation, density, log-density, mean,
  variance, standardized skewness, standardized kurtosis, and raw moments.
- One-component Beta-Skew-t-EGARCH filtering/recursion.
- Two-component Beta-Skew-t-EGARCH filtering/recursion.
- Conditional scale, conditional standard deviation, centered innovations,
  standardized residuals, score values, and latent component paths.
- One- and two-component simulation.
- Exact conditional log-likelihood for both model specifications.
- Bounded maximum-likelihood fitting with configurable asymmetry and skewness.
- Numerical Hessian, inverse observed-information covariance matrix, and
  standard-error extraction.
- BIC per observation, matching the quantity printed by the R package.
- One-step and Monte Carlo multi-step scale and volatility forecasts.
- A simulation/fitting demonstration and a CSV-fitting command-line program.

The regression tests exercise symmetric one-component fitting, full
asymmetric/skewed one-component fitting, and full two-component fitting.

## Deliberate differences from R

- The R package uses `nlminb`; this translation uses a bounded Nelder-Mead
  optimizer with smooth transformations into the constrained parameter space.
  Estimates and convergence diagnostics need not be numerically identical.
- The random-number generator is implemented in Fortran, so a seed does not
  reproduce R's random stream.
- The attached R source appears to contain transcription errors in the central
  third- and fourth-moment formulas. This translation uses the standard central
  moment identities and tests them against simulation.
- In the R multi-step forecast loop, powers are written using the final
  `n.ahead` value for every intermediate horizon. This translation uses the
  current horizon, which is the recursion implied by the documented model.

## Intentionally excluded

Plotting, `zoo` date indexes, package startup behavior, S3 classes, S3 print and
summary dispatch, and other R object infrastructure are not translated. Their
numerical content is available procedurally through the filter and fit result
types.

## Build and test

With GNU Fortran and `make`:

```sh
make
make check
```

`make check` builds with warnings as errors, bounds/runtime checking, and
floating-point traps; runs all numerical tests; executes both applications and
the example; and verifies licensing in every Fortran source file.

An `fpm.toml` manifest is included:

```sh
fpm build
fpm test
fpm run demo_betategarch
```

The validation environment did not contain `fpm`, so those commands are not
claimed as tested. The same sources were compiled directly with GNU Fortran
14.2.0 through the supplied Makefile.

## Programs

Run the demonstration:

```sh
./build/demo_betategarch
```

Fit a CSV column:

```sh
./build/fit_csv FILE [COLUMN=1] [COMPONENTS=1] [ASYM=1] [SKEW=1]
```

For example:

```sh
./build/fit_csv data/example_returns.csv 1 1 0 0
```

Header rows and other rows whose selected field is not numeric are skipped.

## Free-parameter order

One component:

- asymmetric and skewed: `omega, phi1, kappa1, kappastar, df, skew`
- symmetric density: omit `skew`
- no leverage: omit `kappastar`

Two components:

- skewed: `omega, phi1, phi2, kappa1, kappa2, kappastar, df, skew`
- symmetric density: omit `skew`

As in the R package, the two-component model requires asymmetry for
identification.

## Main modules

- `skew_t_mod`: skew-t distribution functions.
- `tegarch_mod`: model types, filtering, simulation, likelihood, fitting,
  covariance calculations, diagnostics, and forecasting.
- `betategarch`: convenience module re-exporting the public API.
