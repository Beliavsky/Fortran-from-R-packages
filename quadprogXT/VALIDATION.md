# Validation

Validation compiler:

```text
GNU Fortran 14.2.0
```

Strict flags:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
```

Regression programs:

1. `test_compact_normalize`
   - reproduces the upstream compact-matrix example;
   - rejects zero constraint columns;
   - verifies unit 2-norm normalization.
2. `test_base_qp`
   - compares `solve_qp_xt` against the supplied base `solve_qp` on the
     canonical 3-variable quadprog example.
3. `test_absolute_value`
   - verifies `sum(abs(b)) <= 1` transformation;
   - verifies an L1 objective penalty against the analytic soft-threshold result.
4. `test_delta_factorized`
   - verifies an `abs(b-b0)` L1 constraint;
   - checks factorized and ordinary expanded QPs agree for `D=I`.
5. `test_full_problem`
   - simultaneously uses ordinary bounds, `abs(b)`, `abs(b-b0)`, and both
     linear penalty types;
   - compares dense and compact solver paths.

All 5 regression programs pass.  Both examples also compile and run.

FPM itself was not installed in the validation environment.  `fpm.toml` was
parsed with Python's TOML parser, and the same source/dependency/test/example
tree was compiled directly with gfortran.
