# Testing

## FPM

```text
fpm build
fpm test
```

The four tests are:

- `test_tail_estimators`: Bessel K reference values, OPP/POP, kurtosis,
  cross-cumulant, Hill, and Pareto estimators;
- `test_elliptical`: Tyler and Cauchy fitting and covariance sanity checks;
- `test_mvt`: Gaussian limit, iterative POP, factor structure, observation
  weights/results, and IEEE-NaN EM imputation;
- `test_mvst`: skewed-t fitting and its reduction to the symmetric Student-t
  model when `gamma = 0`.

## Strict GNU Fortran build

Linux/macOS shell:

```text
./run_gfortran_tests.sh
```

Windows command prompt:

```text
run_gfortran_tests.bat
```

The scripts use Fortran 2018, warnings as errors, runtime bounds/checking, and
backtraces.
