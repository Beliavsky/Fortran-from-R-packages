# Validation

Five deterministic test programs cover:

1. Date arithmetic, holiday rolling, day counts, and the published swap-curve fixture.
2. Positive correlated Black-Scholes factors, deterministic seeding, and rank-2/rank-3 fund mapping.
3. Mortality-factor generation, analytical DBRP/DBRU fixtures, corrected maturity timing, and dispatch through all 19 riders.
4. Equality of one-scenario and duplicated multi-scenario valuation and portfolio aggregation.
5. Synthetic portfolio generation and historical policy aging.

The tests are compiled with:

```text
-std=f2018 -Wall -Wextra -Werror -pedantic -O0 -g
-fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
```

and again with `-O3` and the same warnings-as-errors policy.
