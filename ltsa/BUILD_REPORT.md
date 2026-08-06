# Build report

Compiler: GNU Fortran 14.2.0
Language mode: Fortran 2018

Checked flags:

```text
-O0 -g -std=f2018 -fimplicit-none -Wall -Wextra -Werror
-Wno-compare-reals -fcheck=all -fbacktrace
-ffpe-trap=invalid,zero,overflow -finit-real=snan
```

Optimized flags:

```text
-O3 -march=native -std=f2018 -fimplicit-none
-Wall -Wextra -Werror -Wno-compare-reals
```

Results:

```text
test_durbin: PASS
test_toeplitz: PASS
test_arma_forecast: PASS
test_simulation: PASS
test_innovation: PASS
```

The example completed with finite likelihood, variance, residual, forecast, and
forecast-standard-deviation output. No external BLAS, LAPACK, FFT, R, or C
runtime is required.
