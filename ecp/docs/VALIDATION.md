# Validation

Validation was performed in the translation sandbox with GNU Fortran 14.2.0.

## Direct compiler validation

The maintained sources were compiled in dependency order with strict runtime checking and warnings as errors:

```text
gfortran -std=f2018 -Wall -Wextra -Werror -fcheck=all -fbacktrace
```

Results:

- `test/unit/main.f90`: passed (`All ecp unit tests passed.`)
- `test/algorithms/main.f90`: passed (`All ecp algorithm tests passed.`)
- `test/cp3o_parity/main.f90`: passed (`CP3O pruning parity cases passed, candidates pruned: 6283`)
- `example/change_points.f90`: passed and reported the expected two-change-point path `5 9` for its three-regime fixture.

A second, separate optimized build used:

```text
gfortran -std=f2018 -O2 -Wall -Wextra -Werror
```

All three tests and the example passed again.

The deterministic tests cover all nine exported computational mappings. The divisive permutation branch is seeded explicitly and checked for valid p-values and the requested permutation count without relying on a particular compiler RNG stream. The KCPA tests exercise both default bandwidth estimation and the explicit `bandwidth_rows` interface.

## CP3O optimization regression checks

The release retains `prune=.false.` on all four CP3O entry points as an exhaustive reference path. `test/cp3o_parity` compares pruned and exhaustive results for 12 deterministic multivariate fixtures across all four variants (`e.cp3o`, `e.cp3o_delta`, `ks.cp3o`, and `ks.cp3o_delta`), for 48 comparisons total. Every comparison produced identical GOF vectors, complete change-point path matrices, selected counts, and selected estimates to the stated floating-point tolerance. The fixtures exercised 6,283 actual candidate-pruning decisions.

A separate one-time regression harness compiled the pre-optimization CP3O source from the previous package revision under a renamed module and compared it with the new prefix-sum exhaustive implementation on 10 additional deterministic multivariate fixtures and all four CP3O variants (40 comparisons). GOF vectors and complete path matrices agreed to `1e-12`, confirming that the energy prefix-sum rewrite did not change the exhaustive objective calculation.

For a representative 36-observation, two-variable, four-regime fixture with `K=3`, cached pruning reduced local-statistic evaluations as follows:

| Variant | Pruned | Exhaustive |
| --- | ---: | ---: |
| energy | 867 | 991 |
| delta energy | 846 | 991 |
| KS | 925 | 991 |
| delta KS | 921 | 991 |

These counts are diagnostics rather than timing guarantees. Energy candidates additionally benefit from O(1) prefix-sum evaluation after O(n^2) distance preprocessing.

## KCPA optimization regression check

A separate one-time harness compared the previous direct-block KCPA implementation with the new kernel-prefix implementation on eight deterministic two-variable fixtures. The reconstructed boundary vectors were identical in all cases. The new `segment_costs` implementation computes all contiguous segment costs in O(n^2) after building the kernel matrix, instead of repeatedly summing every kernel submatrix.

## FPM command attempts

The requested literal commands were attempted:

```text
fpm build
fpm test
fpm run --example change_points
fpm clean --all
```

Each returned exit status 127 with `fpm: command not found`, because this sandbox does not provide an FPM executable. These attempts are therefore not reported as successful FPM validation. The manifest itself was parsed successfully as TOML, and its source/test/example layout was exercised by the direct compiler builds above.

## Static audit

The release audit checks:

- exact reconciliation of 9 upstream NAMESPACE exports with 9 translation mappings;
- consistency of totals, translated count, fraction, and empty untranslated list;
- existence and public visibility of every mapped Fortran procedure;
- existence of every recorded upstream and Fortran source path;
- ASCII free-form maintained Fortran and maximum 132-character source lines;
- one dummy argument per declaration, explicit `INTENT` or `VALUE`, and trailing nonempty `!!` FORD documentation;
- exactly one package `dp = real64` definition;
- absence of legacy real kinds or D-exponent literals;
- absence of semicolon-separated maintained Fortran statements;
- absence of self-comparison NaN idioms;
- absence of unsafe fast-math flags and system BLAS/LAPACK links;
- absence of duplicate maintained Fortran source files;
- absence of vendored translated dependencies;
- absence of compiler products, caches, and nested ZIP archives before packaging.

Before archiving, compiler output directories were kept outside the package tree and temporary source backups were removed manually because `fpm clean --all` could not run.
