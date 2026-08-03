# Testing

`fpm test` runs four programs:

- `test_crosscor`: exact regression values from the original C kernels
- `test_cvm`: exact Cramer-von Mises statistics, tables, and combinations
- `test_crossdep`: regression values from `MixedIndTests::EstDepMoebius`
- `test_validation`: invalid dimensions, lags, and degenerate series

For strict GNU Fortran validation:

```text
scripts/test_gfortran.sh
```

Optimized validation can be run with:

```text
FFLAGS="-std=f2018 -Wall -Wextra -Wpedantic -Werror -O3" scripts/test_gfortran.sh
```
