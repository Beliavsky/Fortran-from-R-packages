# NFCP-fortran

A modern Fortran 2018 computational translation of the R package **NFCP**
(N-factor commodity pricing through term-structure estimation).

## Implemented numerical API

- Typed N-factor commodity models with GBM or mean-reverting first factors
- R-compatible parameter layouts and default estimation domains
- `A(T)` and state covariance calculations
- Contract-number and maturity-matched futures stitching
- Dense missing-data Kalman filter, likelihood, filtered states, residuals,
  AIC, and BIC
- Bounded maximum-likelihood estimation using differential evolution followed
  by L-BFGS-B refinement
- Analytical spot and futures forecasts with probability intervals
- Correlated Monte Carlo simulation of spot and futures prices
- Theoretical and empirical futures volatility term structures
- Black-style European options on futures
- Longstaff-Schwartz American/Bermudan options using power, Laguerre, Hermite,
  Legendre, or Chebyshev regression bases

## Build with FPM

```text
fpm build
fpm test
fpm run --example two_factor_oil
```

The bundled `lbfgsb3` package is referenced as a local FPM dependency.

## Minimal model setup

```fortran
use nfcp

type(nfcp_model_t) :: model
call initialize_model(model, n_factors=2, gbm=.true., n_me=1, n_season=0)
model%mu = 0.03_dp
model%mu_rn = 0.01_dp
model%lambda = [0.02_dp, 0.10_dp]
model%kappa = [0.0_dp, 1.2_dp]
model%sigma = [0.20_dp, 0.30_dp]
model%rho = reshape([1.0_dp, -0.2_dp, -0.2_dp, 1.0_dp], [2,2])
model%measurement_error = [0.01_dp]
```

See `docs/API_MAP.md` and `docs/PORTING_NOTES.md` for exact coverage and
intentional differences from R.
