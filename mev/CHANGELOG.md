# Changelog

## 0.3.0

- Add GEV Cox-Snell bias evaluation, implicit correction and Firth-score solve.
- Add the BAB Monte Carlo test and simultaneous envelope.
- Add GPD/GEV tangent-exponential-model profile corrections and modified
  likelihood roots.
- Add uncensored multivariate GP likelihoods for logistic, negative-logistic,
  Brown-Resnick and extremal-t models.
- Add censored multivariate GP likelihoods for the same four model families,
  using native deterministic MV-normal/MV-t QMC for conditional probabilities.
- Add mgp/Poisson/binomial contributions to censored likelihoods.
- Add extremal-t and Brown-Resnick spectral samplers and extend `rmev` to them.
- Add general `rparp` and `rgparp` rejection simulators.
- Correct several clear formula defects in upstream multivariate likelihood
  code; see `PORTING_NOTES.md`.
- Add `test_parity_v03` while retaining all v0.1/v0.2 tests.

## 0.2.0

- Integrate compatible nleqslv and expint Fortran dependency translations.
- Add lower/trimmed Hill and BAB threshold selection.
- Add GPD Cox-Snell/Firth bias correction.
- Add Stein weighted-GPD and exponential-regression tail fitting.
- Add Krupskii-Joe tail-dependence estimation.
- Add deterministic MV-normal/MV-t QMC and multivariate exponent measures.
- Add GPD/GEV profile likelihood routines.
- Add `test_parity_v02` while preserving the v0.1 regression test.

## 0.1.0

- Initial broad numerical translation of mev 2.2.
