# Validation

Validated with GNU Fortran 14.2.0 using:

```text
-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

The permanent test suite contains five programs:

1. `test_reference` -- direct reference case from the original fixed-form NNLS.
2. `test_mixed_sign` -- direct NNNPLS reference with passive-set ordering.
3. `test_exact` -- exact NNLS and mixed-sign recovery.
4. `test_wide_and_rank` -- underdetermined and rank-deficient matrices.
5. `test_status` -- bad-dimension and iteration-limit statuses.

All five pass, along with both examples.

An additional cross-language validation linked the modern modules and the
original `lawson_hanson_nnls.f` / `nnnpls.f` routines in the same executable.
It ran 60 deterministic random matrices of varying `m,n`, once as NNLS and
once with mixed signs (120 solves). Maximum absolute coefficient/residual-norm
difference was approximately `1.1e-14`.

FPM was not installed in the validation environment; the exact FPM source tree
was compiled directly with gfortran, and `fpm.toml` is validated separately.
