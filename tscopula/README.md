# tscopula-fortran

Modern Fortran 2018 translation of the computational code in the R package
`tscopula` 0.3.9. The project is organized as an FPM package and can also be
built with the included Makefile and GNU Fortran scripts.

## Implemented computational scope

- Gaussian ARMA and seasonal ARMA copula processes
- Stationarity and invertibility checks
- Standardized state-space construction, exact stationary Kalman filtering,
  likelihood evaluation, residuals, simulation, prediction, and fitting
- AR/SAR and MA/SMA polynomial expansion
- ACF/PACF transformations and Kendall-PACF generators for ARMA, SARMA(4),
  SARMA(12), ARFIMA, and fractional Brownian noise
- Pair copulas: independence, Gaussian, Student-t, Clayton, Gumbel, Frank,
  Joe, and BB1, including 90/180/270 degree rotations
- Pair-copula CDFs, densities, h-functions, inverse h-functions, Kendall tau
  conversion, and simulation
- Finite-order D-vines, Kendall-PACF-generated D-vines, and D-vines with
  selected non-Gaussian substitutions
- D-vine density and likelihood recursion, Rosenblatt and inverse-Rosenblatt
  functions, conditional prediction, residuals, simulation, and fitting
- Symmetric, degenerate, linear, V2p, V2b, V3p, and V3b V-transforms
- Transform gradients, lower inverses, stochastic inverses, down-branch
  probabilities, and coincidence probabilities
- Optional first-order W-copula branch dependence
- Gaussian, centered Gaussian, Laplace, centered Laplace, skew Laplace,
  double Weibull, skew double Weibull, Student-t, centered Student-t,
  skew Student-t, uniform, and empirical margins
- Parametric margin fitting and Gaussian KDE empirical margins
- Full time-series copula models, IFM-style fitting, likelihoods, simulation,
  residuals, forecasts, quantiles, AIC, BIC, and AICc
- Compatibility interfaces named after the principal upstream routines,
  including `sim`, `fit`, `kendall`, `armacopula`, `sarmacopula`,
  `dvinecopula`, `dvinecopula2`, `dvinecopula3`, `vtscopula`, and `tscm`

## Build

```sh
make test-check
make test-opt
make example
```

With FPM:

```sh
fpm test
fpm run --example tscopula_demo
```

## Minimal example

```fortran
use tscopula

type(arma_copula) :: arma
type(dvine2_copula) :: vine
type(tscm_spec) :: model
real(dp), allocatable :: x(:)

arma = armacopula(ar=[0.55_dp], ma=[-0.20_dp])
vine = arma2dvine('gauss', arma%ar, arma%ma, 4)
model = tscm( &
  vtscopula(tscopula_from_dvine(vine%vine), v2p(0.42_dp,1.25_dp)), &
  margin('sst',[7.0_dp,1.2_dp,0.0_dp,1.0_dp]))
x = sim(model,500)
```

## Design choices

The R package relies on S4 objects, `rvinecopulib`, `FKF`, `ltsa`, `arfima`,
`polynom`, and `kdensity`. The Fortran library replaces those runtime
interfaces with typed derived types and self-contained numerical kernels. No
R, C++, BLAS, LAPACK, or external copula library is required.

See `PORTING_NOTES.md`, `API_MAP.md`, and `VALIDATION.md` for the precise
mapping and numerical-equivalence notes.
