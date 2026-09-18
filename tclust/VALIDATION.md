# Validation

This checkpoint was validated on 2026-09-13 with GNU Fortran 14.2.0.

## FPM availability

FPM is not installed in the validation environment. The required commands were explicitly attempted from the package root and returned exit status 127 (`fpm: command not found`):

```text
fpm build
fpm test
fpm clean --all
```

Consequently, this document does **not** claim that those FPM commands ran successfully on this host. The package includes an FPM manifest, and the same manifest library sources, test program, and example were compiled directly with gfortran as described below.

`fprettify` is also not installed in this environment.

## Strict runtime-check build

All maintained library sources were compiled in module-dependency order together with `test/test_tclust.f90` using:

```text
gfortran -std=f2018 -Wall -Wextra -Werror -pedantic \
  -fcheck=all -fbacktrace -O0 ...
```

Result:

```text
All tclust tests passed.
```

The example compiled and ran with the same flags. The independent parity driver also compiled and ran with the same flags, after which `tools/check_parity.py` reported:

```text
Independent NumPy/SciPy parity checks passed.
```

## Optimized build

The complete deterministic test suite, example, parity driver, and independent Python parity checks were repeated with:

```text
gfortran -std=f2018 -Wall -Wextra -Werror -pedantic -O2 ...
```

All passed.

## Independent numerical checks

`tools/check_parity.py` independently checks major results with NumPy/SciPy rather than reusing the Fortran implementation. It covers:

- adjusted and unadjusted Rand-index calculations;
- Fowlkes-Mallows calculations;
- HARD trimmed-Gaussian classification likelihood from returned centers, covariance matrices, and weights;
- the eigenvalue restriction;
- the trimmed-k-means objective;
- the robust affine-subspace residual objective.

## Static/source audit

`python tools/source_audit.py .` reports:

```text
SOURCE AUDIT PASSED: 16 Fortran files; 11/11 coverage mapping consistent.
```

The audit checks, among other items:

- one package-wide `dp` real kind;
- no legacy `double precision`, `real*8`, D-exponent literals, or `kind(0.0d0)`;
- no executable semicolon-separated statements;
- no self-comparison NaN idioms;
- one dummy argument per declaration with explicit `INTENT`/`VALUE` and a trailing FORD `!!` comment;
- no duplicate Fortran source content;
- README/API/manifest coverage consistency;
- no compiler/build products in the package tree.

The package has no external Fortran dependency and does not copy BLAS, LAPACK, Rcpp, Armadillo, `r.f90`, `r_mod.f90`, or another translated package into the maintained source tree.

## Known compatibility limits

The exported computational-function mapping is 11 of 11, but package status is deliberately **substantial**, not complete. Important differences include the currently untranslated GPCM covariance-pattern branch of `tclust`, R-specific scaling/callback and parallel interfaces, S3/list construction and presentation behavior, a simplified `tclustICsol` solution-ranking heuristic, and non-identical RNG streams.
