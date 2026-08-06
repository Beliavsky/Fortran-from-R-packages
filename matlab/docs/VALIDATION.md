# Validation

The package includes five deterministic test programs:

1. numeric rounding, modular arithmetic, powers, grids, sums, and standard deviations;
2. array construction, flips, rotation, replication, reshape, grids, padding, shape inquiries, and `find`;
3. Hilbert, Vandermonde, magic, Pascal, and Rosser matrices;
4. prime generation, primality, and 64-bit factorization;
5. path parsing/joining, string comparison, and timer operation.

The magic-square tests cover odd order (`n=3,5`), doubly-even order (`n=4`),
and singly-even order (`n=6`), checking row, column, and diagonal sums.

Validated compiler configuration:

```text
GNU Fortran 14.2.0
Fortran 2018
```

Both checked (`-O0 -fcheck=all` plus floating-point traps) and optimized (`-O3`)
builds pass all tests with warnings treated as errors. Exact real comparisons
that intentionally reproduce the upstream source are isolated by
`-Wno-compare-reals`.
