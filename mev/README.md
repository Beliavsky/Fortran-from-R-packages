# mev-fortran

A modern Fortran/FPM translation of the numerical core of the R package
`mev` 2.2 (Modelling of Extreme Values).

Version 0.3.0 closes the main numerical gaps left after v0.2: GEV first-order
bias correction, the BAB Monte Carlo test/envelope, tangent-exponential-model
profile corrections, uncensored and censored multivariate generalized-Pareto
likelihoods, additional Brown-Resnick/extremal-t simulation, and general
R-Pareto/generalized-R-Pareto rejection simulation.

The Fortran API intentionally omits R plotting, S3/formula/model-frame and
presentation infrastructure.

## Implemented numerical areas

- GEV/GPD and seven EGP d/p/q/r families, likelihoods, scores/information,
  fitting, return levels, expected-shortfall and N-block reparameterizations.
- r-largest and Poisson-process threshold likelihoods.
- Hill, lower-trimmed Hill, trimmed Hill, Pickands, moment, de Vries,
  generalized-jackknife, OSZ, generalized-quantile and exponential-regression
  tail-index estimators, second-order rho estimators and Weissman inference.
- BAB threshold selection and the optional Monte Carlo simultaneous-envelope
  test using the translated exponential-integral routines.
- GPD and GEV Cox-Snell/Firth bias correction using the translated nonlinear
  equation solver. GEV Cox-Snell cumulants are evaluated deterministically by
  quadrature of the defining expectation rather than by copying the very large
  generated symbolic expression in the R source.
- Stein weighted-GPD fitting and Krupskii-Joe tail-dependence estimation.
- PWM/L-moments, GPD L-moment fitting, mean residual life and semiparametric
  margins.
- Empirical/Hill tail dependence, extremograms, F-madogram and extremal
  coefficient routines.
- Euclidean/empirical likelihood and empirical Pickands dependence.
- Spatial anisotropy, variograms/correlation and Huesler-Reiss transforms.
- Dirichlet, multivariate normal/t and spectral random generation.
- Logistic, negative-logistic, bilogistic, Dirichlet-mixture,
  Huesler-Reiss, extremal-t and Brown-Resnick spectral samplers.
- Algorithm-1 `rmev` max-stable simulation for the supported spectral models.
- `rparp` and `rgparp` rejection simulation for sum, mean, max, min, L2 and
  site risk functionals (where meaningful for the selected routine).
- Deterministic Halton/antithetic QMC MV-normal and MV-t upper probabilities
  and logistic, negative-logistic, Brown-Resnick, Huesler-Reiss and
  extremal-t exponent measures.
- Uncensored multivariate GP likelihoods for logistic, negative-logistic,
  Brown-Resnick and extremal-t models.
- Censored multivariate GP likelihoods for the same four families, including
  conditional Gaussian/t QMC probabilities and mgp/Poisson/binomial likelihood
  contributions.
- GPD/GEV profile likelihood and TEM modified likelihood-root calculations for
  the principal original parameters.

`API_MAP.md` gives the procedure-level mapping.

## Dependencies and build

Two compatible dependency translations supplied for the v0.2 parity pass are
vendored:

- `nleqslv-fortran` (GPL-2.0-or-later), used by GPD/GEV nonlinear bias
  corrections;
- `expint-fortran` (GPL-3.0-or-later), used by BAB/lower-trimmed-Hill
  threshold-selection code.

`nleqslv-fortran` links BLAS/LAPACK, so the top-level package does too.

```sh
fpm test
fpm run --example demo
```

FPM was not installed in the translation runtime. The FPM source layout was
therefore verified from a clean directory using GNU Fortran 2018. All `mev`
sources compile under `-Wall -Wextra -Werror -fcheck=all`; warnings from the
vendored upstream dependency translations are not promoted to errors.

## License

Upstream `mev` declares `License: GPL-3`; this translation preserves it as
GPL-3.0-only. Compatible vendored dependencies retain their own notices and
licenses under `vendor/`. See `NOTICE.md` and `PORTING_NOTES.md`.
