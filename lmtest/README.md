# lmtest-fortran

A modern free-format Fortran translation of the computational/statistical core of the R package **lmtest 0.9-40**.

The port uses the Fortran Package Manager (FPM), depends only on system BLAS/LAPACK, and keeps the upstream license declaration (`GPL-2 | GPL-3`). The original R source and original fixed-form `pan.f` are retained under `upstream/` for provenance and comparison.

## Implemented numerical functionality

- Ordinary and weighted least squares, fitted values, residuals, covariance matrices, Gaussian log likelihoods, and recursive residuals.
- Coefficient t/z tests and confidence intervals with user-supplied covariance matrices.
- Wald restrictions, nested linear-model F/chi-square comparisons, and likelihood-ratio tests.
- Breusch-Godfrey serial-correlation test, including LM and F versions and the auxiliary-model coefficients/covariance matrix.
- Studentized and nonstudentized Breusch-Pagan heteroskedasticity tests, including observation weights.
- Durbin-Watson statistic with both the upstream Farebrother AS 153 exact-probability algorithm and the lmtest normal approximation.
- Goldfeld-Quandt, Harvey-Collier, Harrison-McCabe, Rainbow, and Ramsey RESET tests.
- Granger-causality testing for two univariate series.
- Cox, Davidson-MacKinnon J, PE, and encompassing tests for nonnested linear models.
- Normal, chi-square, F, and Student-t CDF/survival/quantile support needed by the tests.

## Intentional interface differences

R formula parsing, model frames, S3 methods, automatic model updating, printing/ANOVA formatting, `zoo` indexing, bundled R datasets, and the optional Harrison-McCabe plot are not Fortran numerical algorithms and are not reproduced. Fortran callers pass response vectors and design matrices directly.

For APIs that accept arbitrary R model classes in upstream `lmtest`, the port exposes the underlying numerical operation. For example, `coeftest` becomes `coefficient_tests(beta, vcov, df)`, and `lrtest` becomes `likelihood_ratio_test(loglik1, npar1, loglik2, npar2)`.

## Build

```text
fpm build
fpm test
fpm run --example lmtest_example
```

The manifest links against system `lapack` and `blas`.

## Minimal example

```fortran
use lmtest, only : dp, test_result, breusch_pagan_test
real(dp) :: x(n,p), y(n)
type(test_result) :: bp

bp = breusch_pagan_test(x, y, x)
print *, bp%statistic, bp%p_value
```

See `API_MAPPING.md` for the upstream-to-Fortran mapping and `VALIDATION.md` for numerical validation details.
