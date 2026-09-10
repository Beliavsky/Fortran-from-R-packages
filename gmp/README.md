# gmp - Fortran translation

This directory translates the computational core of the R package **gmp** (upstream version 0.7-5.1) to modern free-form Fortran.  The R package uses the external GNU MP C library; this translation instead contains a package-specific pure-Fortran arbitrary-precision integer/rational engine so that the package has no separately installed GMP, BLAS, or LAPACK requirement.

The retained upstream metadata, R sources, native C/C++ sources, manuals, change log, and MD5 manifest are under `upstream/` for provenance and comparison. They are not compiled by FPM.

## Build

```text
fpm build
fpm test
fpm run --example basic_arithmetic
fpm run --example number_theory
```

The maintained source is standard free-form Fortran 2018 and uses `real(dp)` with `dp = real64` from `gmp_kinds`. No external libraries are required.

## Public API

Use `gmp_api` for the public surface. Major groups include:

- `bigz`: arbitrary-precision signed integers, conversion, exact arithmetic, floor quotient/modulo, GCD/LCM/extended GCD, powers, modular powers/inverses, logs and base-size helpers.
- `bigq`: normalized arbitrary-precision rationals, exact arithmetic, comparison, floor/truncation, ties-to-even rounding, numerator/denominator access.
- number theory and combinatorics: factorial/binomial, Fibonacci/Lucas, Stirling/Eulerian, Bernoulli rationals, exact rational binomial probabilities, primality, next-prime, factorization, deterministic arbitrary-width random integers.
- vectors and matrices: sums/products/cumulative sums/extrema/means, uniqueness/duplication/differences, exact matrix products, transposes, rational solves and inverses.

See `API_COVERAGE.md` and the `[extra.translation]` metadata in `fpm.toml` for limitations and exact R-to-Fortran mappings.

## Important compatibility differences

This is a numerical translation, not an implementation of R's S3/S4 object system or a binary wrapper around GNU MP. In particular, R vector recycling, names/dimnames, replacement/subsetting, NA sentinels, modular-object attributes, formatting/printing, higher-order `apply`, and Rmpfr-dispatched transcendental methods are outside the maintained Fortran API. Primality/factorization and random generation use pure-Fortran algorithms and therefore can differ in performance, probabilistic strategy, or RNG stream from GMP. `solve.bigz` is translated for ordinary integer matrices over the rationals; the upstream modular-matrix solve mode is not.

## Translation coverage

Package status: **substantial**. Coverage is **126 of 132 (95.5%)**. The fraction measures mapped exported computational R functions under the stated coverage basis; it does **not** mean complete R compatibility. Presentation/formatting, data-structure manipulation, R-only generic/interface glue, internal helpers, and functions whose only numerical behavior is delegated to another package are excluded from the denominator.

The six currently unmapped computational functions are the two R NA predicates (`is.na.bigz`, `is.na.bigq`) and the four modular-object attribute getter/setter functions (`modulus`, `modulus.bigz`, `modulus<-`, `modulus<-.bigz`).

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `add.bigz` | `gmp_bigz` | `bigz_add` | substantial |
| `sub.bigz` | `gmp_bigz` | `bigz_sub` | substantial |
| `mul.bigz` | `gmp_bigz` | `bigz_mul` | substantial |
| `div.bigz` | `gmp_api` | `bigz_divide` | substantial |
| `divq.bigz` | `gmp_bigz` | `bigz_quotient` | substantial |
| `%/%.bigz` | `gmp_bigz` | `bigz_quotient` | substantial |
| `mod.bigz` | `gmp_bigz` | `bigz_modulo` | substantial |
| `%%.bigz` | `gmp_bigz` | `bigz_modulo` | substantial |
| `pow.bigz` | `gmp_bigz` | `bigz_pow` | partial |
| `inv.bigz` | `gmp_bigz` | `bigz_invmod` | substantial |
| `!.bigz` | `gmp_api` | `bigz_not` | substantial |
| `|.bigz` | `gmp_api` | `bigz_or` | substantial |
| `&.bigz` | `gmp_api` | `bigz_and` | substantial |
| `xor.bigz` | `gmp_api` | `bigz_xor` | substantial |
| `gcd` | `gmp_bigz` | `bigz_gcd` | substantial |
| `gcd.default` | `gmp_bigz` | `bigz_gcd` | substantial |
| `gcd.bigz` | `gmp_bigz` | `bigz_gcd` | substantial |
| `lcm.default` | `gmp_bigz` | `bigz_lcm` | substantial |
| `lcm.bigz` | `gmp_bigz` | `bigz_lcm` | substantial |
| `as.bigz` | `gmp_bigz` | `bigz_from_int64`, `bigz_from_string` | substantial |
| `as.double.bigz` | `gmp_bigz` | `bigz_to_dp` | substantial |
| `as.integer.bigz` | `gmp_bigz` | `bigz_to_int64` | partial |
| `powm` | `gmp_bigz` | `bigz_powm` | substantial |
| `lt.big` | `gmp_api` | `bigz_lt`, `bigq_lt` | substantial |
| `gt.big` | `gmp_api` | `bigz_gt`, `bigq_gt` | substantial |
| `lte.big` | `gmp_api` | `bigz_lte`, `bigq_lte` | substantial |
| `gte.big` | `gmp_api` | `bigz_gte`, `bigq_gte` | substantial |
| `eq.big` | `gmp_bigz, gmp_bigq` | `bigz_equal`, `bigq_equal` | substantial |
| `neq.big` | `gmp_api` | `bigz_neq`, `bigq_neq` | substantial |
| `is.whole` | `gmp_api` | `real_is_whole`, `bigz_is_whole`, `bigq_is_whole` | substantial |
| `is.whole.default` | `gmp_api` | `real_is_whole` | partial |
| `is.finite.bigz` | `gmp_api` | `bigz_is_finite` | substantial |
| `is.whole.bigz` | `gmp_api` | `bigz_is_whole` | substantial |
| `is.infinite.bigz` | `gmp_api` | `bigz_is_infinite` | substantial |
| `frexpZ` | `gmp_bigz` | `bigz_frexp` | substantial |
| `abs.bigz` | `gmp_bigz` | `bigz_abs` | substantial |
| `sign.bigz` | `gmp_api` | `bigz_signum` | substantial |
| `floor.bigz` | `gmp_api` | `bigz_floor` | substantial |
| `round.bigz` | `gmp_api` | `bigz_round_digits` | substantial |
| `trunc.bigz` | `gmp_api` | `bigz_trunc` | substantial |
| `gamma.bigz` | `gmp_vector` | `bigz_gamma` | partial |
| `cumsum.bigz` | `gmp_vector` | `bigz_cumsum` | substantial |
| `log2.bigz` | `gmp_bigz` | `bigz_log2` | substantial |
| `log.bigz` | `gmp_bigz` | `bigz_log` | partial |
| `log10.bigz` | `gmp_bigz` | `bigz_log10` | substantial |
| `max.bigz` | `gmp_vector` | `bigz_max` | substantial |
| `min.bigz` | `gmp_vector` | `bigz_min` | substantial |
| `prod.bigz` | `gmp_vector` | `bigz_prod` | substantial |
| `sum.bigz` | `gmp_vector` | `bigz_sum` | substantial |
| `duplicated.bigz` | `gmp_vector` | `bigz_duplicated` | substantial |
| `unique.bigz` | `gmp_vector` | `bigz_unique` | substantial |
| `isprime` | `gmp_number_theory` | `isprime_z` | partial |
| `nextprime` | `gmp_number_theory` | `nextprime_z` | substantial |
| `gcdex` | `gmp_bigz` | `bigz_gcdex` | substantial |
| `urand.bigz` | `gmp_number_theory` | `urand_bigz` | partial |
| `sizeinbase` | `gmp_bigz` | `bigz_sizeinbase` | substantial |
| `factorialZ` | `gmp_number_theory` | `factorial_z` | substantial |
| `chooseZ` | `gmp_number_theory` | `choose_z` | substantial |
| `fibnum` | `gmp_number_theory` | `fibnum_z` | substantial |
| `fibnum2` | `gmp_number_theory` | `fibnum2_z` | substantial |
| `lucnum` | `gmp_number_theory` | `lucnum_z` | substantial |
| `lucnum2` | `gmp_number_theory` | `lucnum2_z` | substantial |
| `factorize` | `gmp_number_theory` | `factorize_z` | partial |
| `solve.bigz` | `gmp_matrix` | `bigz_solve`, `bigz_inverse` | partial |
| `add.bigq` | `gmp_bigq` | `bigq_add` | substantial |
| `add.big` | `gmp_bigz, gmp_bigq` | `bigz_add`, `bigq_add` | substantial |
| `sub.bigq` | `gmp_bigq` | `bigq_sub` | substantial |
| `.sub.bigq` | `gmp_bigq` | `bigq_sub`, `bigq_neg` | substantial |
| `sub.big` | `gmp_bigz, gmp_bigq` | `bigz_sub`, `bigq_sub`, `bigz_neg`, `bigq_neg` | substantial |
| `mul.bigq` | `gmp_bigq` | `bigq_mul` | substantial |
| `mul.big` | `gmp_bigz, gmp_bigq` | `bigz_mul`, `bigq_mul` | substantial |
| `div.bigq` | `gmp_bigq` | `bigq_div` | substantial |
| `div.big` | `gmp_api, gmp_bigq` | `bigz_divide`, `bigq_div` | substantial |
| `pow.bigq` | `gmp_bigq` | `bigq_pow` | partial |
| `pow.big` | `gmp_bigz, gmp_bigq` | `bigz_pow`, `bigq_pow` | partial |
| `as.bigq` | `gmp_bigq` | `bigq_make`, `bigq_from_int64`, `bigq_from_string` | substantial |
| `as.double.bigq` | `gmp_bigq` | `bigq_to_dp` | substantial |
| `as.integer.bigq` | `gmp_api` | `bigq_to_int64` | partial |
| `denominator` | `gmp_bigq` | `bigq_denominator` | substantial |
| `denominator<-` | `gmp_api` | `bigq_set_denominator` | substantial |
| `numerator` | `gmp_bigq` | `bigq_numerator` | substantial |
| `numerator<-` | `gmp_api` | `bigq_set_numerator` | substantial |
| `as.bigz.bigq` | `gmp_api` | `bigz_from_bigq` | partial |
| `is.whole.bigq` | `gmp_api` | `bigq_is_whole` | substantial |
| `is.finite.bigq` | `gmp_api` | `bigq_is_finite` | substantial |
| `is.infinite.bigq` | `gmp_api` | `bigq_is_infinite` | substantial |
| `abs.bigq` | `gmp_bigq` | `bigq_abs` | substantial |
| `sign.bigq` | `gmp_api` | `bigq_signum` | substantial |
| `trunc.bigq` | `gmp_bigq` | `bigq_trunc` | substantial |
| `floor.bigq` | `gmp_bigq` | `bigq_floor` | substantial |
| `round0` | `gmp_bigq` | `bigq_round0` | substantial |
| `roundQ` | `gmp_bigq` | `bigq_round_digits` | substantial |
| `round.bigq` | `gmp_bigq` | `bigq_round_digits` | substantial |
| `cumsum.bigq` | `gmp_vector` | `bigq_cumsum` | substantial |
| `mean.bigq` | `gmp_vector` | `bigq_mean` | substantial |
| `solve.bigq` | `gmp_matrix` | `bigq_solve`, `bigq_inverse` | substantial |
| `max.bigq` | `gmp_vector` | `bigq_max` | substantial |
| `min.bigq` | `gmp_vector` | `bigq_min` | substantial |
| `sum.bigq` | `gmp_vector` | `bigq_sum` | substantial |
| `prod.bigq` | `gmp_vector` | `bigq_prod` | substantial |
| `duplicated.bigq` | `gmp_vector` | `bigq_duplicated` | substantial |
| `unique.bigq` | `gmp_vector` | `bigq_unique` | substantial |
| `!.bigq` | `gmp_api` | `bigq_not` | substantial |
| `|.bigq` | `gmp_api` | `bigq_or` | substantial |
| `&.bigq` | `gmp_api` | `bigq_and` | substantial |
| `xor.bigq` | `gmp_api` | `bigq_xor` | substantial |
| `%*%` | `gmp_matrix` | `bigz_matmul`, `bigq_matmul` | substantial |
| `%*%.bigz` | `gmp_matrix` | `bigz_matmul` | substantial |
| `%*%.bigq` | `gmp_matrix` | `bigq_matmul` | substantial |
| `crossprod` | `gmp_matrix` | `bigz_crossprod`, `bigq_crossprod` | partial |
| `crossprod.bigz` | `gmp_matrix` | `bigz_crossprod` | partial |
| `crossprod.bigq` | `gmp_matrix` | `bigq_crossprod` | partial |
| `tcrossprod` | `gmp_matrix` | `bigz_tcrossprod`, `bigq_tcrossprod` | partial |
| `tcrossprod.bigz` | `gmp_matrix` | `bigz_tcrossprod` | partial |
| `tcrossprod.bigq` | `gmp_matrix` | `bigq_tcrossprod` | partial |
| `t.bigz` | `gmp_matrix` | `bigz_transpose` | substantial |
| `t.bigq` | `gmp_matrix` | `bigq_transpose` | substantial |
| `.diff.big` | `gmp_vector` | `bigz_diff`, `bigq_diff` | substantial |
| `Stirling1` | `gmp_number_theory` | `stirling1_z` | substantial |
| `Stirling1.all` | `gmp_number_theory` | `stirling1_all_z` | substantial |
| `Stirling2` | `gmp_number_theory` | `stirling2_z` | substantial |
| `Stirling2.all` | `gmp_number_theory` | `stirling2_all_z` | substantial |
| `Eulerian` | `gmp_number_theory` | `eulerian_z` | substantial |
| `Eulerian.all` | `gmp_number_theory` | `eulerian_all_z` | substantial |
| `BernoulliQ` | `gmp_number_theory` | `bernoulli_q` | substantial |
| `dbinomQ` | `gmp_number_theory` | `dbinom_q` | partial |
