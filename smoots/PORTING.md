# Porting notes

## Source version and scope

This project was translated from `smoots` version 1.1.4. The original R and C++ computational sources are retained in `original/`.

Every exported numerical routine is represented. The following R-specific facilities are intentionally omitted:

- plotting and graphical color helpers
- S3 print, fitted, residual, and plot dispatch
- progress bars
- `future`/`future.apply` parallel process management
- bundled `.rda` demonstration datasets
- R formulas, arbitrary expressions, data frames, and `ts` metadata

## Recovered hidden lookup data

The original package stores essential algorithm tables in `R/sysdata.rda`. They were recovered and translated explicitly:

- the ten `msmooth` algorithm configurations
- trend inflation rates
- derivative inflation rates
- variance-bandwidth multipliers
- fourth-order trend kernels for `p=3`
- first- and second-derivative kernels
- optimal-bandwidth and lower-bound exponents

This avoids replacing the package's data-driven method with generic local regression defaults.

## ARMA estimation substitution

R calls `stats::arima` for:

- `critMatrix`
- AR/MA/ARMA long-run variance estimation
- normal forecasting
- bootstrap re-estimation
- rolling forecast evaluation

A literal translation of `stats::arima` would require porting substantial R `stats` Kalman and optimization infrastructure. The Fortran project instead supplies a self-contained conditional Gaussian estimator:

1. the optional mean is the sample mean;
2. ARMA residuals are generated recursively with zero presample innovations;
3. AR and MA coefficients minimize conditional residual sum of squares;
4. Levenberg-Marquardt steps use an analytical residual Jacobian;
5. Gaussian log-likelihood, AIC, and BIC are calculated from the conditional variance.

This preserves the algorithms and interfaces inside `smoots`, but exact coefficients can differ from R's exact maximum-likelihood estimates. Nonparametric algorithm `A`, `B`, `O`, and `N` paths do not depend on this substitution.

## Numerical integration

The R code approximates kernel integrals using 2,000,001 equally spaced points. The Fortran code uses composite Simpson integration over 20,000 panels. Since every involved kernel is a low-degree polynomial, this is effectively exact to floating-point precision and is considerably cheaper.

## Boundary and invalid-input handling

Several source edge cases were made explicit:

- local-polynomial half-widths are capped so the full window fits the sample;
- singular local regressions return `sm_singular` rather than propagating undefined values;
- lag-window orders are capped at available sample lags;
- confidence-bound square roots are protected against tiny negative roundoff;
- derivative confidence bounds correctly use the pilot trend order for variance estimation;
- `msmooth(method=sm_method_kernel)` forces `p=1`, matching the R wrapper;
- invalid dimensions and probabilities return status codes rather than R exceptions.

## Random numbers and parallelism

R's RNG stream and `future` scheduling are not reproduced. The Fortran bootstrap uses a portable xorshift generator and serial execution. Supplying a seed gives reproducible results across runs of this implementation, but not bit-identical R bootstrap samples.

## `optOrd`

R accepts an arbitrary expression involving row variable `p` and column variable `q`. The Fortran procedure accepts a logical mask matrix. This covers the same computational restriction mechanism without embedding an expression interpreter.
