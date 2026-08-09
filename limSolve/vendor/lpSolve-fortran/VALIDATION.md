# Validation

Compiler used during development:

```text
GNU Fortran 14.2.0
```

Strict flags:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
```

The translated source, all tests, and both examples compile with these flags.
No translated Fortran source/test/example line exceeds the standard 132-column
free-form limit.

## Regression programs

1. `test_general_lp` — upstream continuous LP example: objective 40.5,
   solution `(0, 4.5, 0)`, dual check.
2. `test_mip` — upstream all-integer example: objective 37, solution `(1,4,0)`;
   also verifies three distinct best binary solutions.
3. `test_assignment` — upstream assignment example: objective 8 and exact
   assignment matrix.
4. `test_transport` — Bronson transportation example from upstream docs:
   objective 7790 and expected flow matrix.
5. `test_sparse_q8` — upstream `make.q8` sparse model: 252 nonzeros,
   42 constraints, 64 binary variables, objective 8.
6. `test_statuses` — infeasible and unbounded status detection.

Result: **6/6 pass**.

## Examples

- `basic_lp`
- `assignment`

Both compile and run successfully under the same strict flags.

## FPM

`fpm.toml` parses successfully with Python `tomllib`.  FPM itself was not
installed in the validation container, so the identical FPM source/test/example
tree was compiled directly with gfortran through `scripts/test_gfortran.sh`.

A second strict build/test is performed after extracting the release zip into a
fresh directory, ensuring that no untracked development file is required.
