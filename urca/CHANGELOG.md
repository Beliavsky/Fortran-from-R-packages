# Changelog

## 0.1.0

Initial modern Fortran translation of the computational core of `urca`
1.3-4.

Implemented:

- ADF, ERS, KPSS, Phillips-Perron, Schmidt-Phillips, and Zivot-Andrews tests
- Phillips-Ouliaris cointegration tests
- Johansen trace/maximum-eigenvalue estimation and critical values
- long-run/transitory VECM matrix construction, seasonal and exogenous terms
- `cajools`, `alphaols`, and `cajorls` numerical calculations
- `blrtest`, `alrtest`, `ablrtest`, `bh5lrtest`, `bh6lrtest`, and `lttest`
- `cajolst` endogenous level-shift Johansen procedure
- MacKinnon p-values, quantiles, and unit-root tables
- BLAS/LAPACK regression and linear-algebra support
- FPM manifest, example, and deterministic regression tests

Validation fixes made during translation:

- corrected Zivot-Andrews lagged-difference alignment
- retained upstream-specific `bh6lrtest` iterative beta update
- added Phillips-Perron auxiliary statistics
- guarded MacKinnon response-surface t-ratio evaluation against zero standard
  error without relying on eager/lazy `MERGE` evaluation behavior
- reformatted all source to standard free-form line and continuation limits

The supplied `nlme-fortran` project is retained as a reference but is not a
runtime dependency because the translated `urca` algorithms do not call it.
