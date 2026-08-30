# Validation

Validation was performed in the translation environment on 2026-08-29.

## Dependency/API validation

The maintained package imports:

- `dp` and `i64` from `r_kinds` in sibling `rfortran-core`.
- `r_median` from `r_quantiles` in sibling `rfortran-core`.
- `r_mad` from `r_robust` in sibling `rfortran-core`.
- `symmetric_eigen` from `r_linalg` in sibling `rfortran-linalg`.

The current repository interfaces were checked directly. `r_kinds::dp` is `real64`; the `r_linalg` symmetric eigensolver also accepts `real(real64)`, so the shared linear-algebra boundary uses the same real kind as this package's `dp`.

## Strict compiler validation

Because the execution sandbox does not contain FPM, the source was additionally compiled with GNU Fortran 14.2.0 using validation-only interface-compatible stubs located outside the release tree for the shared modules, and system LAPACK/BLAS for the MDS eigensolver.

Compiler options:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

Results:

```text
all randomForest tests passed
training errors: 0
first-case class probabilities: 1.0000 0.0000
training MSE: 0.053346
final OOB MSE: 0.352393
```

The test suite covers numeric and categorical classification/regression, stratified sampling, OOB/permutation importance, proximities, unsupervised forests, margins/outliers/class centers, rough imputation, proximity-based iterative imputation, `mtry` tuning, `rfcv`, partial dependence, and classical MDS distance reconstruction.

`fpm.toml` was parsed as TOML and its sibling dependency paths were verified.

## Release hygiene

The pre-archive scanner reported:

```text
maintained Fortran files: 12
unique source hashes: 12
maximum maintained Fortran line length: 128
extra style issues: []
```

It checks for copied shared/dependency source, duplicate Fortran source, build/cache/archive artifacts, semicolon-separated statements, forbidden real-kind forms, D-exponent literals, alternate `real64`/`iso_fortran_env` kinds in maintained code, and lines longer than 132 characters.

## Required FPM command attempts

Immediately before archive creation, the requested commands were run literally from the package root:

```text
$ fpm build
bash: fpm: command not found
exit status 127

$ fpm test
bash: fpm: command not found
exit status 127

$ fpm clean --all
bash: fpm: command not found
exit status 127
```

Therefore this environment cannot truthfully claim a successful FPM invocation. No `build/` directory or FPM artifact was created. `tools/release_check.py` performs the exact `fpm build`, `fpm test`, hygiene scan, and `fpm clean --all` sequence on an FPM-equipped checkout.
