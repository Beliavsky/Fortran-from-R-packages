# Testing

The test suite covers every one of the 55 exported calculations using the examples supplied in the original `.Rd` files, plus:

- compatibility-generic dispatch;
- vector-valued residual-income routines;
- mean and median selector branches;
- GGM/comparable and EBITDA/sales selector branches;
- round-half-to-even behavior;
- array-size mismatch handling with IEEE NaN; and
- independent unrounded valuation identities.

Recommended commands:

```text
fpm test
fpm run
fpm run --example reference_examples
```

For a strict direct build with GNU Fortran:

```text
gfortran -std=f2018 -Wall -Wextra -Werror -pedantic \
  -fcheck=all -ffpe-trap=invalid,zero,overflow -O0 \
  src/stock_analyst.f90 test/test_stock_analyst.f90 -o test_stock_analyst
./test_stock_analyst
```
