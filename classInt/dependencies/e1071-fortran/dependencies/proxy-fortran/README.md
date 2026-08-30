# proxy-fortran

Modern Fortran/FPM translation of the computational core of the R package `proxy` 0.4-29.

The library provides auto-, cross-, and row-paired proximity calculations; the built-in numeric, binary, nominal, mixed/Gower, Mahalanobis, and Levenshtein measures; an extensible typed registry for custom Fortran proximity callbacks; and helpers for R-compatible packed `dist` storage.

## Basic use

```fortran
use proxy, only: dp, proxy_dist_auto

real(dp), allocatable :: d(:, :)
real(dp) :: x(3, 2)

x(1, :) = [0.0_dp, 0.0_dp]
x(2, :) = [1.0_dp, 0.0_dp]
x(3, :) = [1.0_dp, 1.0_dp]
call proxy_dist_auto(x, "Euclidean", d)
```

The public `proxy_dist_auto`, `proxy_dist_cross`, and `proxy_dist_pairwise` routines return distances for numeric matrices. If the named method is natively a similarity, the package applies the same conversion convention as R `proxy`. Corresponding `proxy_simil_*` routines convert native distances to similarities.

String edit distances use `levenshtein_distance`, `levenshtein_auto`, and `levenshtein_cross` because Fortran character arrays and numeric matrices are statically distinct types.

## Built-in methods

Translated built-ins include Euclidean/L2, Mahalanobis, Bhjattacharyya, Manhattan/L1, supremum/Chebyshev, Minkowski/Lp, Canberra, Wave/Hedges, divergence, Kullback-Leibler, Bray-Curtis, Soergel, Podani, chord, geodesic, Whittaker, Hellinger, fuzzy Jaccard, cosine, angular, extended Jaccard, extended Dice, correlation, all binary coefficients registered by upstream `proxy`, nominal chi-square based coefficients, Gower similarity, and Levenshtein edit distance.

`mutual_information_similarity` is also retained from the upstream C computational source even though current upstream registration no longer exposes that C routine through the package registry.

The translation intentionally preserves a few non-obvious upstream formula semantics. In particular, the R implementations of Wave/Hedges and Soergel use scalar `min()`/`max()` calls across their vector arguments rather than componentwise `pmin()`/`pmax()`, and the matrix fast path for extended Dice uses the factor of two present in the upstream C routine.

## Gower data

For mixed variables, encode rows in a real matrix and supply one variable-type code per column to `gower_auto_similarity` or `gower_cross_similarity`:

- `proxy_gower_logical`: zero/one logical data
- `proxy_gower_factor`: numeric factor-level codes
- `proxy_gower_metric`: ordinary numeric variables
- `proxy_gower_ordinal`: ordinal integer codes transformed with the same internal-code convention as upstream

NaN represents missing numeric data. Optional ranges, minima, and weights correspond to the computational arguments used by upstream Gower preprocessing.

## Extensible registry

Custom numeric and binary proximity functions can be added with `proxy_register_numeric` and `proxy_register_binary`. The callback interfaces are public. Fortran procedure pointers are used instead of R's dynamic function objects.

## Packed `dist` utilities

`proxy_pack_dist` and `proxy_unpack_dist` use the same strict-lower-triangle order as R's `dist` class. `proxy_subset_dist`, `proxy_row_sums_dist`, `proxy_row_means_dist`, `proxy_row_dist`, and `proxy_col_dist` provide the corresponding numerical utilities without R attributes or character subscripting.

## Build

With FPM:

```text
fpm test
fpm run --example basic
```

The release is also validated directly with GNU Fortran using strict standard, warning, bounds, interface, and floating-point checks.

All maintained reals use the single public `dp` kind from `proxy_kinds`. Every maintained dummy argument has explicit `INTENT`/`VALUE`, is declared on its own declaration line, and has a meaningful trailing FORD `!!` comment. Maintained source is free-form and `fprettify`-compatible.
