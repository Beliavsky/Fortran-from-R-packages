# Validation

The release was compiled with gfortran 14.2 using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O2
```

All six permanent test programs pass:

1. `test_criteria` -- distance, MST, phi-p, and all seven discrepancy criteria.
2. `test_classic` -- factorial, LHS, OLH, NOLH/NOLHDR, Faure, scale/unscale.
3. `test_stochastic` -- D-max, Strauss, and WSP invariants.
4. `test_lhs_opt` -- SA/ESE criterion improvement and Latin-hypercube invariants.
5. `test_uniformity` -- uniformity statistics/quantiles and 2-D/3-D RSS values.
6. `test_nolh_all` -- every NOLH/NOLHDR dimension from 2 through 29 and Latin
   level structure.

The example `basic_designs` also compiles and runs under the same flags. In the
release build it reduced the C2 discrepancy of a centered 12 x 3 LHS from
0.086169 to 0.075026.

## Independent/randomized checks

- 200 deterministic pseudo-random designs were evaluated independently in
  NumPy for all seven discrepancy measures plus minimum distance, coverage,
  mesh ratio, and phi-p. Maximum absolute difference: approximately
  1.8e-15.
- The translated 2-D and 3-D projection-CDF kernels were compared against the
  original `src/C_rss.c` on 10,000 identical numeric input quadruples, including
  deliberately near-degenerate coefficients around the upstream 1e-12 cutoff.
  All 20,000 CDF values were bit/numerically identical in the printed double
  representation (maximum absolute difference 0).
- Fixed deterministic RSS regression cases are included for both 2-D and 3-D
  designs, including their worst projection directions.

`fpm` itself was not installed in the validation environment, so the exact FPM
source/test/example tree was compiled directly with gfortran. `fpm.toml` is
included and separately parsed as TOML during release validation.
