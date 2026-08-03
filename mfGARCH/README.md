# mfGARCH-fortran

A self-contained modern Fortran translation of the computational core of the
R package `mfGARCH` 0.2.2. It estimates, forecasts, and simulates multiplicative
GARCH-MIDAS models in which

```text
conditional variance = long-run MIDAS component * short-run GARCH component.
```

The project uses Fortran Package Manager conventions and has no external
runtime dependencies.

## Implemented functionality

- Restricted and unrestricted normalized beta-MIDAS lag weights.
- One or two low-frequency explanatory variables.
- Irregular mixed-frequency mappings through integer period-index arrays.
- GARCH(1,1) and asymmetric GJR-GARCH(1,1) short-run components.
- Gaussian log likelihood and observation-level likelihood contributions.
- Constrained estimation through smooth parameter transformations.
- BFGS and Nelder-Mead optimization with optional multi-stage refinement.
- Numerical Hessian, OPG covariance, and sandwich/robust covariance.
- Fitted `tau`, `g`, standardized residuals, BIC, variance ratio, and next-period
  long-run variance forecast.
- Multi-horizon conditional-variance forecasts.
- Standard GARCH-MIDAS simulation with Gaussian or standardized Student-t
  intraday innovations.
- Realized-variance-dependent GARCH-MIDAS simulation.
- Andersen-style diffusion-limit simulation.
- Reusable low-level kernels corresponding to the upstream C++ routines.

Plotting, R S3 methods, data-frame manipulation, and packaged R data objects
are intentionally excluded.

## Build and test

With FPM:

```sh
fpm build
fpm test
fpm run demo_mfgarch
```

Without FPM, GNU Fortran build scripts are supplied:

```sh
./scripts/test_gfortran.sh
./scripts/test_gfortran_optimized.sh
```

On Windows with `gfortran` available in `PATH`:

```bat
scripts\test_gfortran.bat
```

## Mixed-frequency input

The high-frequency return vector and an integer period-index vector have the
same length. Each period index points to the corresponding element of the
low-frequency covariate vector. For example, five daily observations per week
can be represented by

```fortran
period = [1,1,1,1,1, 2,2,2,2,2, ...]
```

This representation also supports irregular periods because the number of
high-frequency observations assigned to each low-frequency period need not be
constant.

## Minimal fitting example

```fortran
use mfgarch

type(mfgarch_model) :: start
type(mfgarch_fit_result) :: fit
integer :: status

start%k = 12
start%asymmetric = .true.
start%alpha = 0.02_dp
start%beta = 0.85_dp
start%gamma = 0.04_dp
start%m = 0.0_dp
start%theta = 0.0_dp
start%w1 = 1.0_dp
start%w2 = 3.0_dp

call fit_mfgarch(returns, period, start, fit, status, covariate=macro_series)
call print_fit_summary(fit)
```

For numerical stability, the upstream package recommends scaling log returns
by 100 before fitting. The same advice applies here.

## References

- Engle, R. F., Ghysels, E., and Sohn, B. (2013). Stock market volatility and
  macroeconomic fundamentals. *Review of Economics and Statistics*, 95(3),
  776-797.
- Conrad, C. and Kleen, O. (2020). Two are better than one: Volatility
  forecasting using multiplicative component GARCH models. *Journal of
  Applied Econometrics*, 35(1), 19-45.
