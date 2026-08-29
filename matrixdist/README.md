# matrixdist-fortran

Modern Fortran/FPM translation of the computational core of the R package **matrixdist 1.1.9** by Martin Bladt and Jorge Yslas, with Alaric Mueller as contributor.

The port uses numerical arrays and derived types instead of reproducing R's S4 object system.  The umbrella module is `matrixdist`.

## Implemented numerical areas

- Continuous PH: density/CDF/survival/hazard/quantile, moments, Laplace/MGF, likelihood, simulation, sums/mixtures/minima/maxima, random structures, and censored EM fitting.
- Inhomogeneous PH: Weibull, Pareto, matrix-lognormal, log-logistic, Gompertz and GEV transformations; density/CDF/simulation; fixed-transform likelihoods; scale handling.
- Discrete PH: mass/CDF/survival/PGF, factorial moments, mean/variance, likelihood, simulation, operations and EM fitting.
- Feed-forward bivariate PH/DPH: density/tail, transforms, moments/covariances, simulation, likelihood and EM fitting.
- Conditional multivariate PH/DPH: joint density/CDF/PGF/Laplace, mixed moments, mean/covariance, simulation, likelihood, multivariate DPH EM and right-censored multivariate PH EM.
- Inhomogeneous bivariate/multivariate transforms.
- MPH*: reward-process mean/covariance, simulation, time-value-of-reward transformations, linear combinations, reward sanitization, and the reward-matrix EM update.
- Survival-regression numerical evaluation for `reg`, `reg2`, and `aft` parameterizations.
- Matrix exponential, Van Loan blocks, integer powers, inverses/solves, Kronecker products/sums, block diagonals, and uniformization helpers.
- LRT and covariance-to-correlation helpers.

The supplied MIT-licensed `r_mod.f90` is used for the common R-compatible RNG/probability helper layer rather than duplicating it.  The exact original is retained in `upstream/r_mod-original.f90`; the build copy differs only by free-form line continuation formatting inherited from the previously validated helper copy.

## Deliberately omitted R-facing layers

Formula/model-frame construction, S4 classes/method dispatch, plotting, printing, progress output, `reshape2` data reshaping, and presentation-only methods are omitted.  The high-level `nnet::multinom` mixture-of-experts orchestration and generic `stats::optim` outer loops for estimating transformation/regression coefficients are not recreated as R-like object APIs; the PH/DPH EM, likelihood, regression-evaluation, and conditional-responsibility mathematics they operate on are available as array-oriented kernels.

See `API_MAPPING.md` and `PORTING_NOTES.md` for the exact translation boundary and documented upstream edge corrections.

## Build

```text
fpm build
fpm test
fpm run --example basic_matrixdist
```

Linear systems and matrix inverses use the repository's shared
`rfortran-linalg` dependency and its pinned pure-Fortran LAPACK backend.
System BLAS and LAPACK libraries are not required.

## License

`matrixdist`-derived code is distributed under GPL-3.0-only, matching the upstream package declaration.  `r_mod.f90` remains MIT licensed.  See `NOTICE.md` and `licenses/`.
