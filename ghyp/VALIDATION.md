# Validation

## Toolchain

The release was compiled with GNU Fortran 14.2.0 in both checked and optimized
configurations.

Checked flags:

```text
-std=f2018 -O0 -Wall -Wextra -Wconversion-extra -Werror
-fcheck=all -fbacktrace
```

Optimized flags:

```text
-std=f2018 -O2 -Wall -Wextra -Wconversion-extra -Werror
```

## Test suites

```text
test_gig: PASS
test_distributions: PASS
test_risk_transform: PASS
test_fitting: PASS
test_portfolio: PASS
```

The demo and both examples also compile and run in checked and optimized builds.

## Independent numerical references

The tests include independent SciPy and NumPy reference values for:

- real-order modified Bessel K;
- GIG density, CDF, quantile, mean, and variance;
- generalized-hyperbolic density, CDF, and quantiles;
- generalized-hyperbolic mean, covariance, skewness, and Pearson kurtosis;
- multivariate generalized-hyperbolic density;
- expected shortfall and Omega ratio;
- Gaussian and Student special cases;
- analytical minimum-variance and target-return portfolios.

The principal GH reference uses

```text
lambda = 0.7
alpha  = 1.8
delta  = 1.2
beta   = 0.3
mu     = 0.2
```

and agrees with SciPy's `genhyperbolic` implementation to the tolerances encoded
in the tests.

## Portability

Assertions use combined scale-aware tolerances rather than machine-epsilon-only
comparisons. This avoids false failures from minor differences in `erfc`,
`log_gamma`, and floating-point expression evaluation across platforms.

The implementation does not require nested-procedure trampolines and links
without GNU executable-stack warnings.

## Release audit

The release audit checks:

- valid TOML syntax;
- ASCII-only translated text;
- free-form source lines no longer than 132 characters;
- `implicit none` in every Fortran compilation unit;
- SPDX licensing headers in every translated Fortran file;
- original and translated SHA-256 manifests;
- a clean rebuild and test run from the exact ZIP archive.
