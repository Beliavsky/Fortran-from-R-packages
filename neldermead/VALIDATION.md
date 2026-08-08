# Validation

Validation was performed with GNU Fortran 14.2.0 using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
```

The strict build compiles all library source files without warning suppression,
then builds and runs every test and example.

## Tests

Six test programs pass:

1. `test_fminsearch` -- variable-shape Nelder-Mead on Rosenbrock.
2. `test_fixed` -- Spendley-Hext-Himsworth fixed-shape search.
3. `test_box` -- bound-constrained Box complex method.
4. `test_constraints` -- nonlinear constrained Michalewicz G6 problem; with the
   deterministic regression seed the solution is approximately
   `(14.0950001, 0.8429609)` with objective `-6961.8137`.
5. `test_grid` -- bounded Cartesian grid search and result ordering.
6. `test_simplex` -- Spendley simplex construction and evaluation accounting.

Both examples also build and run successfully.

The objective, constraint, output, and termination callbacks are declared through
explicit abstract interfaces. Calls do not rely on implicit external interfaces.

On GNU/Linux, linking tests/examples that pass *internal* procedures as callback
actual arguments may emit an executable-stack trampoline note. That is a linker
implementation detail for contained-procedure callbacks and is not an implicit-
interface diagnostic from the translated library.

FPM was not installed in the validation environment. `fpm.toml` was parsed with
a TOML parser, and the same `src/`, `test/`, and `example/` tree was compiled and
run directly with gfortran.
