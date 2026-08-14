# Validation

Validated with gfortran using:

```text
-std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

Permanent tests cover:

1. the upstream-style "choose the k smallest values" problem;
2. subset range, cardinality, uniqueness, and final row sorting;
3. monotonic best-history behavior under elitism;
4. exact `mutfrac` to `mutprob` conversion;
5. supplied initial populations and summary calculations; and
6. 20 independently generated `n=12, k=4` additive problems whose exact global
   optima are known by sorting all component costs.

All 20 randomized exact-optimum cases converge to the analytical optimum with
the fixed regression seeds/settings in `test_randomized_exact.f90`.

The example also converges to subset `[1,2,3]` with objective `0.3`.

The FPM manifest is parsed during the release audit. If `fpm` is unavailable,
the exact source/test/example tree can be compiled directly with the flags
above.
