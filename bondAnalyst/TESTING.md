# Testing

The test suite contains two programs.

## `test_all_exports`

Covers every one of the 55 original exports using the package's documented
examples. Expected values were calculated independently from the published
formulas. It also tests the two added conventional helpers and status handling.

## `test_identities`

Checks:

- YTM repricing after six-decimal output rounding;
- Z-spread repricing after four-decimal output rounding;
- Macaulay/modified-duration identities;
- the settlement-period full-price example;
- duration calculated from an externally supplied full price;
- spot/forward compounding consistency;
- invalid day counts, impossible positive-root cases, and size mismatches.

## Compiler configurations used

Strict checked build:

```text
gfortran -std=f2018 -Wall -Wextra -Werror -pedantic \
  -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
```

Optimized build:

```text
gfortran -std=f2018 -O3 -Wall -Wextra -Werror -pedantic
```

Both test programs, the application, and the example are rebuilt from clean
source trees in both configurations.
