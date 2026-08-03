# Testing

## FPM

```text
fpm test
```

## Direct GNU Fortran validation

```text
scripts/test_gfortran.sh
scripts/test_gfortran_optimized.sh
```

The strict script uses:

```text
-std=f2018 -Wall -Wextra -Wpedantic -Werror
-fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
```

The optimized script uses `-O3` and keeps all warnings as errors.

## Test programs

### `test_gaussian`

Checks exact low-penalty coefficient recovery, coefficient bounds, exclusions,
dense/CSC equivalence, predictions, coefficient interpolation, nonzero sets,
cross-validation, and relaxed refitting.

### `test_glm_families`

Checks binomial signs and AUC, Poisson fit and deviance, multiresponse Gaussian
coefficients, grouped multinomial classification, the matrix-response
multinomial interface, assessment, and multinomial cross-validation.

### `test_cox`

Checks Cox path fitting, concordance, analytical partial-likelihood gradients
against central finite differences, Efron tied events, strata, and Cox
cross-validation.

### `test_utilities`

Checks the custom-family callback, missing-value replacement, constant-column
removal, matrix concatenation, deterministic multinomial sampling, ROC/AUC,
confusion tables, measure names, survival-data construction, and error paths.
