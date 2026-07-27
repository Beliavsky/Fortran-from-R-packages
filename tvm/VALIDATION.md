# Validation

## Environment

The release was validated with GNU Fortran 14.2.0 using Fortran 2018 mode.
FPM was not installed in the validation container, so the same source graph was
compiled directly in FPM module dependency order. `fpm.toml` was parsed and the
standard `src`, `app`, `example`, and `test` layouts were audited.

## Strict build flags

```text
-std=f2018
-Wall
-Wextra
-Wpedantic
-Wconversion-extra
-Wimplicit-interface
-Werror
-fcheck=all
-fbacktrace
-O0
```

A second build was run with `-O2`.

## Test programs

```text
test_cashflows: PASS
test_curves: PASS
test_interpolation: PASS
test_rate_curve: PASS
```

The translated project contains 1527 lines of Fortran across the library, tests, demo, and examples.
The demo and both examples also compile and run.

## Numerical coverage

Validation includes:

- fixed reference values for NPV, irregular NPV, IRR, XIRR, payment, total
  financial cost, remaining balance, and discount-spread adjustment;
- bullet, French, German, and grace-period loan schedules;
- upstream one-rate, flat-curve, and upward-curve identities;
- futures, effective-zero, nominal-zero, continuously compounded, swap,
  French-loan, and German-loan transformations;
- direct/effective and discount/rate round trips;
- monotone interpolation and exact node recovery;
- rate scaling, callback curves, and curve-based present value.

## Linker note

GNU Fortran on Linux may report that `tvm_cashflows.o` requests an executable
stack. This results from passing internal objective procedures to the generic
root solver. It does not affect test results. Compilers that implement internal
procedure callbacks without trampolines do not emit this warning.

## Reproduction

Run the portable FPM commands:

```text
fpm build
fpm test
fpm run tvm_demo
fpm run --example basic_cashflows
fpm run --example rate_curves
```

Or run `scripts/validate.sh` on a Unix-like system with GNU Fortran.
