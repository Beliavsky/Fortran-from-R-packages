# Computational coverage

## Included

| Upstream area | Fortran implementation |
|---|---|
| `dgig`, `pgig`, `qgig`, `rgig` | Complete density, CDF, quantile, and random generation |
| `Egig`, `ESgig` | Raw moments, mean, variance, log mean, inverse mean, and ES |
| `ghyp`, `ghyp.ad` | Typed chi/psi and alpha/delta constructors |
| `hyp`, `NIG`, `student.t`, `VG`, `gauss` | Specialized constructors |
| `dghyp`, `pghyp`, `qghyp`, `rghyp` | Univariate distribution API and multivariate density/simulation |
| Multivariate `pghyp` | Simulation-based rectangle probability with standard error |
| `ghyp.moment` | Exact integer moments and numerical absolute/noninteger moments |
| `ghyp.skewness`, `ghyp.kurtosis` | Standardized central moments |
| `ESghyp`, `ghyp.omega` | Tail means and Omega ratio |
| `ESghyp.attribution` | Finite-difference ES sensitivities and Euler contributions |
| `transform`, `scale`, `[` | Linear transformation, standardization, and marginal extraction |
| `coef(..., alpha.delta)` | Alpha/delta conversion through `ghyp_alpha_delta` |
| `qqghyp` calculations | Plot-independent QQ arrays through `qqghyp_data` |
| `fit.ghypuv` and family wrappers | Native direct maximum likelihood for all six families |
| `fit.ghypmv` and family wrappers | Native direct multivariate maximum likelihood |
| `stepAIC.ghyp` | AIC comparison of GH, hyp, NIG, Student, VG, and Gaussian |
| `lik.ratio.test` | Likelihood-ratio statistic and chi-square p-value |
| `portfolio.optimize` | SD, VaR, and ES optimization for minimum-risk, tangency, and target-return problems |
| `mean`, `vcov`, distribution naming | Typed moment and family-name routines |

## Replaced implementation details

- R's adaptive `integrate` calls are replaced by cached Gauss-Legendre
  quadrature and explicit infinite-domain transformations.
- R's `besselK` is replaced by a native real-order Bessel K integral.
- The compiled upstream GIG rejection sampler is replaced by direct gamma or
  inverse-gamma generation in limiting cases and inverse-CDF generation in the
  general case.
- The specialized multivariate EM fitter is replaced by a common direct
  likelihood optimizer.
- R's `optim` and `numDeriv` dependencies are replaced by native Nelder-Mead
  and numerical Hessians.

## Excluded

The following are presentation or R-runtime infrastructure rather than
standalone numerical algorithms:

- S4 class registration and dispatch;
- plotting, histogram, line, pairs, and attribution charts;
- formatted show/summary methods;
- R formula, list, and data-frame adapters;
- compiled use of bundled example data;
- R vignette generation and documentation machinery.
