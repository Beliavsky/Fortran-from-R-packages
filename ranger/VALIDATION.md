# Validation record

This file records the release checks performed for the modern Fortran translation of `ranger` 0.18.0.

## Toolchain available in the translation sandbox

- GNU Fortran: `GNU Fortran (Debian 14.2.0-19) 14.2.0`
- Fortran standard used for the independent validation build: Fortran 2018
- FPM: not installed in this sandbox

The package has one sibling FPM dependency, `rfortran-core`. The maintained `ranger` source only imports `r_kinds::dp` and `r_kinds::i64` from that dependency. The independent build used the current `r_kinds` interface (`dp = real64`, `i64 = int64`) and did not place a copy in the package tree.

## Strict independent compiler validation

A fresh external build directory was used. Every maintained library source was compiled with:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fimplicit-none -fcheck=all -O0 -g
```

The full deterministic test program then compiled, linked, and ran successfully:

```text
all ranger tests passed
```

Both examples also compiled, linked, and ran successfully. Representative output was:

```text
training errors: 0
first case vote fractions:    1.00000   0.00000
survival for first case:    0.9273   0.8009   0.5165   0.5165   0.5165   0.5165
```

The tests exercise classification, probability, regression, quantile-regression preparation/prediction, survival forests, categorical and missing-value handling, sampling/holdout behavior, multiple ranger split rules, OOB behavior, importance calculations, hierarchical shrinkage, infinitesimal-jackknife variance, and forest utilities. They include class-specific zero node/bucket limits, survival `time_interest` normalization/count selection, `oob_error=.false.` with permutation importance still enabled, and the zero-OOB-coverage NaN behavior of ranger.

This checkpoint also adds strict-runtime regression coverage for Maxstat cutpoint counts/unadjusted tails and the IJ rule that calibration is disabled for 20 or fewer prediction points. The strict build exposed and fixed a survival factor-order insertion-sort bounds defect and several expressions that relied on C/C++-style short-circuit assumptions; the maintained Fortran now branches before potentially out-of-range array references. Additional parity work accepts zero entries in class-specific node/bucket vectors, reproduces ranger's scalar and vector `time.interest` preprocessing, implements the `oob.error` computation switch in all four forest engines, and aligns zero-OOB prediction/error sentinels.

## Source and release hygiene

The pre-archive scan verified:

- 19 maintained free-form `.f90` files under `src/`, `test/`, and `example/`;
- 19 distinct SHA-256 hashes for those maintained files;
- maximum maintained Fortran line length of 132 characters;
- all maintained Fortran files import the shared `dp` from `r_kinds`;
- no `double precision`, `real*8`, `kind(0.0d0)`, direct `real64`, package-local `iso_fortran_env`, or D/d-exponent real literals in maintained code;
- all maintained decimal/exponent real constants have `_dp` kind suffixes;
- no semicolon-separated Fortran statements;
- no copied `rfortran-core`, BLAS, LAPACK, ARPACK, Rcpp, RcppEigen, or Matrix source;
- no duplicate Fortran source files;
- no object files, module files, executables, libraries, caches, build directories, or nested ZIP archives in the release tree;
- `fpm.toml` parses as TOML.

`API_COVERAGE.md` records implementation-level adaptations rather than silently claiming bit-for-bit equivalence. Corrected impurity importance and ranger's Lau92/Lau94 Maxstat corrections are translated; the principal remaining numerical adaptation is that empirical-Bayes IJ calibration uses a native two-parameter Newton/line-search optimizer instead of R's `nlm`.

## Required FPM gate

The required commands were attempted literally in the package root before packaging. This sandbox has no `fpm` executable, so they cannot be reported as successful:

```text
$ fpm build
bash: line 1: fpm: command not found
exit_status=127

$ fpm test
bash: line 1: fpm: command not found
exit_status=127

$ fpm clean --all
bash: line 1: fpm: command not found
exit_status=127
```

`tools/release_check.py` runs `fpm build`, `fpm test`, the maintained-source/style/provenance hygiene checks, and finally `fpm clean --all` on an FPM-equipped checkout. The absence of FPM here is an environment limitation, not a claimed FPM validation pass.

## Post-archive check

After the ZIP is created, it is extracted into a fresh external directory and the strict GNU Fortran build, deterministic test suite, examples, and archive-layout/hygiene checks are repeated. The final ZIP SHA-256 is reported alongside the delivered artifact.
