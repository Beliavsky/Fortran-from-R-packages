# Validation

## Compiler configuration

The project was validated with GNU Fortran 14.2.0 in two configurations.

Checked build:

```text
-std=f2018 -Wall -Wextra -Werror -Wconversion-extra
-Wimplicit-interface -fcheck=all -fbacktrace -O0
```

Optimized build:

```text
-std=f2018 -Wall -Wextra -Werror -Wconversion-extra
-Wimplicit-interface -O2
```

## Test programs

- `test_distribution`: fixed SciPy references, CDF/quantile inversion, density
  integration, flags, macro factors, boundaries, and seeded simulation moments.
- `test_fit`: fixed NumPy/SciPy OLS and transformed-parameter references,
  residual identities, bias correction, finite-portfolio correction, and errors.
- `test_prediction`: fixed conditional link, mean, median, and tail-quantile
  references for new macroeconomic observations.
- `test_inference`: fixed IID delta-method covariance, exact symmetry, native HAC
  covariance, and confidence-interval width checks.

Expected output:

```text
test_distribution: PASS
test_fit: PASS
test_inference: PASS
test_prediction: PASS
```

The demo and both examples are also compiled and executed.

## Independent references

Distribution values and deterministic fit, prediction, and covariance
references were calculated independently with NumPy and SciPy. The density was
also integrated numerically to one.

## FPM

FPM was not installed in the validation environment. The manifest parses as
TOML and the project follows FPM's standard automatic `src`, `test`, `app`, and
`example` target discovery. Direct validation scripts are included for Linux
and Windows.
