# Validation

Compiler used during translation:

```text
GNU Fortran 14.2.0
```

## Checked build

The checked validation script uses:

```text
-std=f2018
-O0
-Wall
-Wextra
-Wconversion-extra
-Wimplicit-interface
-Werror
-fcheck=all
-fbacktrace
```

All tests passed:

```text
test_copula_risk: PASS
test_distributions: PASS
test_fitting: PASS
test_mixtures: PASS
test_skewt: PASS
validation: PASS
```

## Optimized build

The same tests passed at `-O2` with warnings promoted to errors:

```text
optimized validation: PASS
```

## Reference coverage

- Multivariate normal and Student densities were checked against independent
  SciPy values.
- Student-copula density, 99% VaR, 99% expected shortfall and tail dependence
  were checked against independent SciPy formulas.
- Normal, Student and gamma quantile/CDF inversions were checked numerically.
- Constant and inverse-gamma Mahalanobis-mixture identities were checked against
  chi-square and F laws.
- Student and grouped-Student simulations were checked against theoretical
  variances and grouped correlation formulas.
- Zero-skew skew-t density/CDF were checked against Student formulas.
- Skew-t simulation was checked against its theoretical mean.
- Normal and Student fitting were tested on deterministic synthetic samples.
- QQ output was checked for valid sorted observed distances.

## Reproduction

Linux/macOS:

```text
./scripts/validate.sh
./scripts/validate_optimized.sh
```

Windows with GNU Fortran:

```text
scripts\validate.bat
```

FPM was not available in the translation environment. `fpm.toml` uses standard
automatic discovery for `src`, `app`, `example` and `test` targets.
