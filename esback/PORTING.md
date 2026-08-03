# Porting notes

## esreg dependency

`esback::esr_backtest` delegates joint VaR/ES regression and covariance estimation to `esreg`. The Fortran project embeds the numerical subset used by esback rather than requiring a separate R runtime.

Translated elements include:

- Fissler-Ziegel loss for `g1=2`, `g2=1`
- strict, auxiliary, and intercept ESR design matrices
- quantile-regression starting values
- Nelder-Mead local search with deterministic multistarts
- density-at-quantile estimation (`iid` and `nid`)
- conditional location-scale maximum likelihood
- conditional truncated variance (`ind`, `scl_N`, `scl_sp`)
- lambda, sigma, and sandwich covariance matrices

## Numerical substitutions

- R's `quantreg::rq` is replaced by direct minimization of the check loss.
- R's `optim` is replaced by a self-contained Nelder-Mead implementation.
- R's `density(..., bw="SJ")` path is represented by a stable Gaussian plug-in bandwidth and the same trapezoidal truncated-moment calculation. This may cause small covariance and p-value differences from R for `sigma_scl_sp`.
- R's random-number generator is replaced by an explicit portable xorshift generator. Bootstrap results are deterministic for a given seed but are not expected to reproduce R's exact sample stream.
- Dense inverses and solves use LAPACK/BLAS.

## Preserved source behavior

- `er_backtest` restarts its bootstrap stream from seed 1 for the simple and standardized tests.
- Hommel and Bonferroni corrections follow the R formulas exactly.
- The executable `esreg` source's `g2=1` second derivative uses `-1/z^3`; this source convention is retained even though its documentation states `-2/z^3`.
- The default ESR covariance is misspecification robust.

## Omitted infrastructure

Formula parsing, S3 methods, data frames, printing, plotting, and the serialized `risk_forecasts.rda` loader are omitted. Callers pass numeric vectors/matrices directly.
