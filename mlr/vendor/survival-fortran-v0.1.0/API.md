# API

All floating-point APIs use `real(dp)` from `survival_kinds`.

## Nonparametric survival

- `kaplan_meier(time, status, fit [, weights])`
- `kaplan_meier_counting(start, stop, status, fit [, weights])`
- `aalen_johansen(time, from_state, to_state, nstate, prob, utime [, initial])`

`survfit_result` contains `time`, `n_risk`, `n_event`, `n_censor`, `survival`,
`cumhaz`, `std_err`, and `std_chaz`.

## Cox proportional hazards

- `coxph_fit(time, status, x, result [, method, weights, offset, maxiter, eps])`
- `coxph_fit_counting(start, stop, status, x, result, ...)`
- `cox_baseline(time, status, x, beta, base_time, cumhaz [, method, weights])`
- `cox_martingale_residuals(...)`
- `cox_schoenfeld_residuals(...)`

`method` is `"breslow"` or `"efron"`.

## Parametric AFT

- `survreg_fit(time, status, x, dist, result [, weights, maxiter, eps])`
- `survreg_loglik(...)`

Implemented names: `extreme`, `weibull`, `exponential`, `rayleigh`, `gaussian`,
`lognormal`, `loggaussian`, `logistic`, `loglogistic`.

## Tests/statistics

- `survdiff(time, status, group, ngroup, result [, rho, strata])`
- `concordance_right(time, status, risk, result [, weights])`

## Data transformations

- `finegray_expand(tstart, tstop, ctime, cprob, extend, keep, ...)`
- `surv_split(start, stop, status, cut, ...)`
- `pseudo_survival(time, status, eval_time, pseudo)`

## Penalized splines

- `pspline_basis(x, nterm, degree, boundary, intercept, basis, penalty, knots [, status])`

This uses the attached `splines-fortran` dependency and reproduces survival's
equally-spaced knot construction plus the second-difference penalty matrix.
