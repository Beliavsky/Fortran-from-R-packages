# Testing

The package includes four test programs:

- `test_tvm`: future/present values, annuities, NPV, IRR, discount-rate solving, payments, periods, and perpetuities
- `test_rates_ratios`: yield conversions and financial ratios
- `test_accounting`: EPS, diluted EPS, weighted shares, issuable shares, depreciation, and FIFO/LIFO/WAC inventory costing
- `test_statistics`: means, sampling error, weighted portfolio return, and TWRR

Run with FPM:

```console
fpm test
```

The supplied compiler-only script builds with strict diagnostics and runtime checks:

```console
./run_gfortran_tests.sh
```

The default GNU Fortran flags are:

```console
-std=f2018 -Wall -Wextra -Werror -pedantic \
-fcheck=all -ffpe-trap=invalid,zero,overflow \
-ffree-line-length-none -O0 -g
```

Applications and examples are also compiled and executed by the script.
