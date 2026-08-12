# Validation

Validated with GNU Fortran 14.2.0.

## Strict debug build

The library and all tests compile with:

```text
-std=f2018 -O0 -g -fcheck=all
-Wall -Wextra -Werror
-Wimplicit-interface -Werror=implicit-interface
```

All five FPM regression programs pass:

```text
test_objective_io: PASS
test_swap_delta: PASS
test_sa_had20: PASS
test_reproducibility: PASS
test_qaplib_formats: PASS
```

## Optimized build

The same suite passes at `-O2`.

## QAPLIB sweep

A separate validation driver called `read_qaplib` on every bundled `.dat`
file. All 136 files parsed successfully.

The reader intentionally does not require that every raw `.sln` permutation
reproduce its recorded optimum under one fixed permutation orientation because
some historical QAPLIB files use the inverse assignment convention. This is
also the behavior of the original R reader.

## Numerical checks

`had12`:

```text
known optimum = 1652
qap_obj(solution) = 1652
```

`had20`, 10 restarts, Fortran seed 1000:

```text
known optimum = 6922
best objective = 6922
```

Every pairwise swap in a six-variable deterministic regression is also checked
by comparing `qap_swap_delta` with a complete O(n^2) objective recomputation.
