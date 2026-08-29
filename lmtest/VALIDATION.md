# Validation

Validated with GNU Fortran 14.2.0 using:

```text
-std=f2018 -Wall -Wextra -Werror=implicit-interface -fcheck=all
```

and system LAPACK/BLAS.

The regression suite covers:

1. Published Farebrother AS 153 probability benchmarks from the comments in upstream `pan.f`.
2. Normal, Student-t, chi-square, and F probabilities/quantiles against independent SciPy reference values.
3. OLS coefficients, covariance-derived standard errors, log likelihoods, coefficient tests, and confidence intervals.
4. Breusch-Godfrey LM/F, Breusch-Pagan (both forms), Goldfeld-Quandt, Harvey-Collier, Rainbow, RESET, and Durbin-Watson on a deterministic regression dataset. Reference values were independently calculated with Statsmodels/SciPy or direct implementations of the upstream formulas.
5. Cox, J, PE, and encompassing tests on deterministic nonnested model pairs.
6. Granger causality on a deterministic bivariate recursive series.

All included test programs pass under runtime bounds/checking.

FPM was not installed in the validation container, so the source graph was compiled directly with `gfortran`; `fpm.toml` is included for normal FPM builds and links to `lapack`/`blas`.
