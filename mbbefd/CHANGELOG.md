# Changelog

## 0.1.0

- Initial modern Fortran/FPM computational port of mbbefd 0.8.14.
- Port both MBBEFD parameterizations and associated exposure/moment functions.
- Port shifted truncated Pareto and generalized beta first-kind laws.
- Port generic and named one-inflated distributions.
- Port empirical total-loss/EECF statistics and Theil helpers.
- Port MLE and TLM moment fitting plus parametric/nonparametric bootstrap.
- Use supplied `fitdistrplus-fortran`, `actuar`, and `alabama` ports as FPM dependencies.
- Omit plotting, S3 display methods, and Rcpp glue.
- Correct generic one-inflated quantile handling and two upstream `(a,b)` constraint/Jacobian inconsistencies.
