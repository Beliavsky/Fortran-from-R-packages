# API coverage: TMB 1.9.25 → Fortran

## Implemented

| Upstream area | Fortran API | Status |
|---|---|---|
| `distributions_R.hpp` normal helpers | `dnorm`, `pnorm`, `qnorm` | Implemented scalar API |
| exponential helpers | `dexp`, `pexp`, `qexp` | Implemented |
| Weibull helpers | `dweibull`, `pweibull`, `qweibull` | Implemented |
| binomial mass | `dbinom`, `dbinom_robust` | Implemented |
| beta/gamma/lognormal/logistic density | `dbeta`, `dgamma`, `dlnorm`, `dlogis` | Implemented |
| F and Student-t density | `df_density`, `dt_density` | Implemented |
| skew-normal density | `dsn` | Implemented |
| multinomial mass | `dmultinom` | Implemented |
| sinh-asinh family | `dshasho`, `pshasho`, `qshasho`, `norm2shasho` | Implemented |
| `tmbutils/kronecker.hpp` | `kronecker_product` | Implemented dense real matrices |
| `tmbutils/order.hpp` | `order_real`, `sort_real` | Implemented deterministic ascending order |
| interval lookup used by spline utilities | `find_interval` | Implemented |
| `tmbutils/interpol.hpp` | `interpolate2d` | Implemented value interpolation; NaNs omitted |
| `tmbutils/density.hpp` MVN | `mvnorm_nll` | Implemented dense SPD covariance |
| `tmbutils/density.hpp` N01 | `n01_nll` | Implemented |
| `tmbutils/density.hpp` AR1 | `ar1_nll`, `ar1_mvn_nll` | Implemented scalar and dense multivariate Gaussian forms |
| `tmbutils/density.hpp` UNSTRUCTURED_CORR | `unstructured_corr`, `unstructured_corr_nll` | Implemented |
| `tmbutils/matexp.hpp` | `matrix_exponential` | Implemented dense scaling/squaring Taylor form |
| `tmbutils/romberg.hpp` | `romberg_integrate` | Implemented one-dimensional Romberg integration |
| derivative facility | `gradient_fd`, `hessian_fd` | Numerical finite differences; not exact AD |

## Not yet translated / not claimed as parity

- TMBad and CppAD tape construction, reverse/forward automatic differentiation, atomic operators,
  graph transforms, code generation, and sparsity discovery.
- `MakeADFun`, `sdreport`, R `.Call` wrappers, dynamic compilation/loading, RStudio/debug helpers,
  plotting/profile interfaces, and other R-specific workflow code.
- Eigen expression-template API and CHOLMOD sparse factorization/inverse-subset machinery.
- Higher structured-density classes (`GMRF`, `SEPARABLE`, `PROJ`, `ARk`, continuous AR2) beyond the currently exposed dense MVN, unstructured-correlation, and AR1 kernels.
- Specialized distributions/functions requiring larger numerical support (incomplete gamma/beta
  derivatives and inverses, Poisson CDF, Bessel families, Tweedie, COM-Poisson, COM-binomial).
- Spline objects and derivative-aware atomic spline operators.
- TMB simulation macros and random-number integration with R.

These are meaningful remaining parity targets; the present package is a portable computational
subset rather than a replacement for TMB's C++ AD framework.
