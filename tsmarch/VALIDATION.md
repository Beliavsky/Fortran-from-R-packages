# Validation report

## Compiler

- GNU Fortran 14.2
- Fortran 2018 free-form source

## Checked configuration

The checked build uses:

- `-std=f2018`
- `-Wall -Wextra -Werror -pedantic`
- `-fcheck=all -fbacktrace`

All five test programs pass:

1. `test_dcc`: DCC simulation, dynamic filtering, constant correlation,
   normalization, and positive conditional variances
2. `test_copula`: Gaussian and Student transforms and copula filtering
3. `test_ica_gogarch`: whitening, FastICA/RADICAL reconstruction, GO-GARCH
   covariance/correlation, higher moments, and portfolio moments
4. `test_utils_risk_fft`: EWMA/Ledoit-Wolf covariance, multivariate densities
   and draws, combinations, ESCC, VaR/ES, and convolution distribution methods
5. `test_estimation_workflow`: end-to-end DCC and GO-GARCH estimation and
   simulation-based forecasting

## Optimized configuration

The optimized build uses `-O3` with strict warning and standards checks.  The
same five programs pass without compiler warnings.

## Demonstration

`example/demo_tsmarch.f90` simulates three correlated standardized innovation
series, scales them to returns, estimates a constant-correlation multivariate
model with EWMA marginals, and generates a five-step Monte Carlo forecast.

The validated run produced finite likelihoods, positive volatility forecasts,
and a positive-definite fitted correlation matrix.

## Reproducibility limits

- Random streams differ from R.
- ICA signs and ordering are non-unique.
- Finite-difference optimization and inference need not duplicate R's solver
  path exactly.
- Direct DFT convolution is deterministic but not identical to every FFT
  normalization/rounding path.
