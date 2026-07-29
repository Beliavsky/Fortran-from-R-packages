# BEKKs-fortran

A modern Fortran/FPM translation of the computational algorithms in the R package **BEKKs 1.4.7**, for multivariate BEKK conditional-volatility modelling.

The project provides a typed API for estimation, simulation, filtering, forecasting, volatility impulse responses, value at risk, forecast backtesting, Monte Carlo evaluation, and the package's matrix-algebra utilities. The original R and C++ package tree is retained under `original/` for license and provenance auditing.

## Implemented models

- Full BEKK(1,1)
- Diagonal BEKK(1,1)
- Scalar BEKK(1,1)
- Asymmetric versions of all three models
- User-defined sign vectors for the asymmetric indicator

For a return vector `r(t)`, the full symmetric recursion is

```text
H(t) = C C' + A' r(t-1) r(t-1)' A + G' H(t-1) G
```

The asymmetric model adds

```text
I(t-1) B' r(t-1) r(t-1)' B
```

where `I(t-1)` is one when every component satisfies the selected sign condition.

## Main numerical features

- Gaussian quasi-maximum likelihood
- Berndt-Hall-Hall-Hausman optimization
- Per-observation score matrices
- Numerical likelihood Hessians
- OPG and QML sandwich covariance estimates
- Deterministic and randomized starting-value searches
- Stationarity testing using the spectral radius of the Kronecker transition
- Conditional covariance filtering and standardized residuals
- Fixed-innovation and random simulation
- Multi-step covariance, standard-deviation, and correlation forecasts
- Parameter-uncertainty forecast bands
- Volatility impulse-response functions with numerical delta-method bands
- Marginal and portfolio VaR using normal, empirical, or standardized Student-t residual quantiles
- Kupiec and Christoffersen coverage tests
- Rolling-window VaR backtests
- Multivariate portmanteau testing
- Monte Carlo parameter-recovery evaluation
- Elimination, duplication, commutation, selection, and cut matrices
- `vech`, inverse-`vech`, lag-matrix, covariance extraction, symmetric square roots, and generalized inverses

## Data layout

Time-series matrices use shape `T x N`:

- rows are observations;
- columns are assets or series.

Conditional covariance paths use shape `N x N x T`.

## Build with FPM

BLAS and LAPACK are required.

```bash
fpm build
fpm test
fpm run
fpm run --example asymmetric_bekk
fpm run --example portfolio_var
```

## Reproducible GNU Fortran build

When FPM is unavailable:

```bash
./scripts/test_gfortran.sh all
```

The script performs clean debug and optimized builds, runs every test, and compiles and runs every application and example.

## Minimal example

```fortran
use iso_fortran_env, only: int64
use bekks
implicit none

type(rng_state) :: rng
type(bekk_spec_type) :: spec
type(bekk_fit_result) :: fit
real(dp), allocatable :: data(:,:), h(:,:,:), theta(:)
integer :: status

! theta must contain lower-triangular C followed by model coefficients.
call rng_seed(rng,12345_int64)
call simulate_sbekk(theta,1000,2,rng,data,h,status)
spec=bekk_spec(bekk_scalar,.false.,initial_theta=theta)
call bekk_fit(spec,data,fit,max_iter=50,use_qml=.true.)
```

See `app/`, `example/`, and `API.md` for complete working programs and procedure signatures.

## Scope

Plotting, S3 methods, `xts`/`zoo`/`ts` metadata, progress displays, and R parallel orchestration are not translated. Their underlying numerical data are returned through typed Fortran result objects. The original bundled R datasets are retained under `original/data/` but are not parsed by the Fortran library.

Important numerical and behavioral differences are described in `PORTING.md`.

## License and citation

The original package and this derivative translation use the MIT License. See `LICENSE` and `NOTICE.md`.

The original package citation is:

M. J. Fuelle, A. Lange, C. M. Hafner, and H. Herwartz (2024), "BEKKs: An R Package for Estimation of Conditional Volatility of Multivariate Time Series," *Journal of Statistical Software*, 111(4), 1-34.
