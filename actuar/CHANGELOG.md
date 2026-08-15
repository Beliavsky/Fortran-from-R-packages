# Changelog

## 0.3.0

- Added minimum-distance estimation for individual/grouped CvM, modified
  chi-square and layer-average-severity criteria.
- Added deductible/franchise/limit/coinsurance/inflation coverage
  transformations for CDFs and mixed density/mass functions.
- Added the Hachemeister barycenter regression-credibility algorithm with
  weighted orthogonalization and original-basis coefficient recovery.
- Ported the exact iterative hierarchical-credibility recursion from
  `src/hierarc.c`, with Buhlmann-Gisler and Ohlsson alternatives.
- Added a callback-driven hierarchical portfolio simulator corresponding to
  the numerical role of `rcomphierarc`.
- Added `test_v03` and `v03_remaining` example.
- Retained the upstream R/C sources used for the v0.3 translation under
  `upstream/reference-v03`.

## 0.2.0

- Added Feller-Pareto, inverse Pareto and inverse transformed-gamma families.
- Added broad heavy-tail limited-moment coverage and upstream incomplete-beta
  continuation.
- Vendored the previously translated `expint-fortran` dependency for extended
  incomplete-gamma calculations.
- Corrected the inverse-exponential scale parameterization from v0.1.
- Added exact aggregate convolution, normal/normal-power approximations and
  compound/mixture simulation helpers.
- Added phase-type Cramer-Lundberg and general Sparre-Andersen ruin algorithms.
- Added hierarchical credibility and Hachemeister regression credibility.
- Added grouped moments, quantiles, ogive and empirical limited expected value
  calculations.
- Added beta/chi-square and inverse-distribution supplementary moments.
- Added `test_v02` and `v02_extended` example.

## 0.1.0

- Initial modern Fortran/FPM translation of the computational core of
  `actuar` 3.3-7.
- Added heavy-tailed loss distributions and actuarial moment helpers.
- Added zero-truncated and zero-modified frequency laws.
- Ported the upstream Poisson-inverse-Gaussian probability recurrence.
- Added phase-type density, CDF, moments, MGF and simulation.
- Added Panjer recursion, discretization, aggregate VaR/CTE and ruin bound.
- Added Buhlmann-Straub and selected Bayesian credibility calculations.
- Added strict GNU Fortran regression tests and example program.
