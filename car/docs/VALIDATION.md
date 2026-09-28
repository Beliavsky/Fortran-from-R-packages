# Validation

Validation performed in the translation environment used GNU Fortran 14.2.0.

The full library plus deterministic test program was compiled from a clean external build directory with:

```text
gfortran -std=f2018 -Wall -Wextra -pedantic -fcheck=all -fbacktrace \
  src/car_kinds.f90 src/car_utils.f90 src/car_transforms.f90 \
  src/car_contrasts.f90 src/car_geometry.f90 src/car_kernels.f90 \
  src/car_matrix.f90 src/car_diagnostics.f90 src/car_data.f90 \
  src/car_power_fit.f90 src/car_box_tidwell.f90 src/car_po_test.f90 \
  src/car_lm.f90 src/car_anova.f90 src/car_manova.f90 \
  src/car_hypothesis.f90 src/car.f90 test/test_car.f90 -o test_car
./test_car
```

Result: `test_car: PASS`.

The example was also compiled and executed against the same source set. The deterministic suite includes the Box-Tidwell and proportional-odds diagnostics, linear-model summary/confidence/prediction/bootstrap helpers, Type III LM ANOVA, multivariate Type III SSP and MANOVA statistics, coefficient matching, linear-hypothesis parsing, exact name lookup, deterministic sampling-index injection, character-to-factor coding, and coefficient-comparison statistics. See the recorded validation notes in the release/archive creation log if one accompanies this tree.


Additional release checks verified that:

- `fpm.toml` parses as TOML and contains 51 unique function mappings for a 51-function coverage denominator;
- `README.md`, `docs/API_COVERAGE.md`, and `fpm.toml` list the same 51 mapped R functions;
- every maintained Fortran dummy data argument has explicit `INTENT` or `VALUE` and a trailing `!!` FORD comment;
- no dummy procedure arguments remain, avoiding the Fortran restriction that prevents `INTENT` on a procedure dummy;
- maintained Fortran source contains no statement-separating semicolons, legacy double-precision kind idioms, or non-ASCII text; and
- the packaged source tree contains no compiler products, caches, nested ZIP archives, or other build artifacts.

The sandbox used to create this translation did not have the Fortran Package Manager executable installed, and outbound package installation was unavailable. Therefore an actual `fpm build` / `fpm test` invocation could not be performed in that environment. The manifest was parsed as TOML and the exact source/test/example sets were compiled directly with gfortran. The intended downstream verification remains:

```text
fpm build
fpm test
fpm run --example basic
```
