# Rfast2-fortran

Modern Fortran/FPM computational-core translation of **Rfast2 0.1.5.6**.

This project uses the supplied **Rfast-fortran 0.3.0** translation as an FPM
path dependency. It does not require R, Rcpp, RcppArmadillo, or RcppParallel.

## Implemented Rfast2-specific areas

- Type-7 scalar/row/column quantiles and trimmed means
- intersection/merge and matrix-triangular helpers
- column grouping summaries
- jackknife means and column means/variances
- binary and continuous prediction metrics
- portable PCG32 uniforms, sampling, beta/gamma/exponential/chi-square,
  geometric, Cauchy, t and related RNG helpers
- Jarque-Bera tests, Pearson correlation tests, covariance helpers,
  empirical entropy, pooled variances
- circular correlations, PINAR(1), Kaplan-Meier, Moran's I
- Poisson rate Wald tests, Walter intervals, WLS/REML meta-analysis,
  silhouette, permutation/bootstrap t tests and univariate energy tests
- Rfast2-specific univariate MLEs including gamma-Poisson, half-Cauchy,
  zero-centered Cauchy, Kumaraswamy, power-law, zero-inflated gamma,
  zero-inflated Laplace/Weibull, simplex, generalized-normal-zero,
  unit-Weibull, continuous Bernoulli and SP
- column-wise beta/Cauchy/half-Cauchy/lognormal/logit-normal/Borel/power-law/
  unit-Weibull/SP/half-normal fitting
- constrained, robust and heteroscedastic least squares
- logistic, Poisson, gamma, Weibull, multinomial, grouped-binomial,
  zero-truncated Poisson and tobit regression plus univariate scan routines
- PCA, PCR, Mahalanobis depth, leverage, item difficulty/discrimination,
  and covariance-distance calculations

The vendored Rfast dependency additionally supplies its extensive reusable
array, linear-algebra, distance, distribution, MLE, regression, testing,
directional, naive-Bayes, score-test, repeated-measures, variable-selection,
and PC-skeleton numerical kernels.

## Build

```bash
fpm build
fpm test
```

or compile with a Fortran 2018 compiler. The project itself has no C/C++
source dependency.

## Example

```fortran
program demo
   use rfast2
   implicit none
   real(dp) :: x(5)

   x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
   print *, quantile_rfast2(x, 0.25_dp)
end program demo
```

## Scope

This is a computational-core v0.1 port, not a claim that every one of
Rfast2's roughly 199 exported R entry points has a one-for-one Fortran
wrapper. `docs/TRANSLATION_COVERAGE.md` lists the high-value remaining
numerical targets. R formula/S3/data-frame/environment/parallel orchestration
layers are intentionally out of scope.
