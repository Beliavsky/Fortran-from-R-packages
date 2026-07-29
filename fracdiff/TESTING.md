# Testing

## FPM

```sh
fpm test
fpm run
fpm run --example long_memory_estimators
fpm run --example custom_innovations
```

## GNU Fortran audit

```sh
./scripts/test_gfortran.sh
```

The script performs two complete builds:

1. strict Fortran 2018 with warnings as errors, bounds/runtime checks,
   backtraces, and floating-point traps;
2. optimized `-O3` compilation.

Every test, application, and example is rebuilt in both configurations.

## Test coverage

### Fractional differencing

- independent fixed convolution reference
- FFT/direct agreement
- integer first-difference identity
- fractional-weight recursion
- polynomial multiplication and AR root calculation

### Haslett-Raftery filtering and ARMA derivatives

- fixed filtered-series reference
- estimated generalized mean
- innovation log-variance sum
- analytical ARMA Jacobian against central finite differences

### Simulation

- fixed-innovation source recursion
- post-filter mean handling
- deterministic seeding
- finite and plausible generated values

### Semiparametric estimators

- independent fixed GPH estimate and both standard errors
- independent fixed Sperio estimate and both standard errors

### Maximum likelihood and inference

- ARFIMA(1,d,1) parameter recovery
- fractional-noise parameter recovery
- innovation-scale recovery
- symmetric Hessian and covariance
- positive standard errors
- confidence-interval construction
- AIC/BIC consistency
- covariance recomputation with a changed finite-difference interval
