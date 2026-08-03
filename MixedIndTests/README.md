# MixedIndTests-fortran

Modern Fortran translation of the computational code in the R package
`MixedIndTests` 1.2.0. The library implements tests of independence and
randomness for continuous, discrete, or mixed data using multilinear empirical
copulas, multiplier bootstrap statistics, pairwise rank dependence, and
Moebius copula-based covariances.

The project is self-contained and uses only standard Fortran 2018. It can be
built with FPM:

```text
fpm build
fpm test
fpm run
```

The public module is `mixedindtests`. The original R and C sources are retained
under `original/` and the translated source is under `src/`.

## Included computations

- empirical unique values, CDF, and probability masses;
- multilinear Kendall and Spearman dependence coefficients;
- pairwise and serial portmanteau tests;
- nonserial, univariate-serial, and vector-serial Cramer-von Mises kernels;
- Gaussian multiplier bootstrap and Fisher combinations;
- Spearman, van der Waerden, and Savage Moebius coefficients;
- data-driven lag-order selection;
- Poisson AR(1) and copula-series simulation;
- seven marginal inverse-distribution functions used by the upstream package.

Plot-only functions `AutoDep`, `Dependogram`, and `DependogramZ` are omitted.
Parallel R worker management is replaced by deterministic sequential bootstrap
loops. See `PORTING.md` for numerical and interface details.

## Validation

The tests include exact regression values generated from the original C source
for scalar and vector serial statistics, multiplier matrices, and Moebius
coefficients. Run `scripts/test_gfortran.sh` or the Windows batch equivalent for
a strict compiler build independent of FPM.
