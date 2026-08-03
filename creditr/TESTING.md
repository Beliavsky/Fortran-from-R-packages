# Testing

## Configurations

The release is tested in two clean configurations:

### Strict/debug

```text
-std=f2018 -Wall -Wextra -Werror -pedantic
-fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
```

### Optimized

```text
-std=f2018 -O3 -Wall -Wextra -Werror -pedantic
```

Run either configuration with:

```text
./scripts/build_gfortran.sh debug
./scripts/build_gfortran.sh release
```

## Test programs

- `test_dates_and_formulas`: date conventions, holidays, simple probability/spread formulas, implied recovery, and PV01.
- `test_curve`: CSV input and fixed ISDA zero-curve nodes.
- `test_cds_reference`: four fixed CDS cases against the retained ISDA C engine.
- `test_conversion_and_risk`: spread/upfront inversion, spread DV01, recovery risk, and CS10.

Every application and example is also compiled and executed by the script.

## Reference tolerances

Money-market curve points are checked at machine precision. Swap-bootstrapped points use absolute tolerances of `2e-12` to `3e-12`. CDS cash values and sensitivities use relative tolerances appropriate for independent Fortran and C implementations of iterative integrations.
