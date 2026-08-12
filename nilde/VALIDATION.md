# Validation

Translation environment: GNU Fortran 14.2.0.

## Strict build

The library, all five tests and both examples compile with:

```text
-std=f2018 -O0 -g -fcheck=all
-Wall -Wextra -Werror
-Wimplicit-interface -Werror=implicit-interface
```

Passing tests:

1. `test_nlde`
   - 10 solutions for `3*s1 + 2*s2 + 5*s3 + 16*s4 = 18` with the default part bound.
   - 3 solutions when exactly six parts are required.
   - binary equality regression.
   - 108 binary inequality solutions for the package documentation example.
2. `test_partitions_subsetsum`
   - 20 partitions of 8 into at most six parts.
   - 2 partitions of 8 into exactly six positive parts.
   - documented 0/1 and bounded subset-sum examples.
3. `test_knapsack`
   - bounded and 0/1 optimum/tie behavior.
   - original unbounded-return compatibility and corrected maximizing mode.
4. `test_binpacking`
   - weights `(70,60,50,40,30,20,10)`, capacity 100.
   - minimum 3 bins and 9 canonical optimal packings.
5. `test_tsp`
   - assignment lower bound 4.
   - exact Hamiltonian optimum 22.
   - 19 integer objective levels tested and 4 optimal directed tours.

## Optimized build

The complete five-test suite also passes with:

```text
-std=f2018 -O2 -Wall -Wextra -Werror
-Wimplicit-interface -Werror=implicit-interface
```

## Source format

No Fortran source line exceeds the standard 132-column free-form limit.
The FPM manifest parses as TOML. FPM itself was not installed in the translation
environment, so equivalent direct GNU Fortran builds were used.
